"""Build a portable fixture; retain separate tests for legacy Speedict banks."""

from pathlib import Path

import pytest

from seqbank import SeqBank


@pytest.fixture(scope="session")
def portable_bank(tmp_path_factory):
    path = tmp_path_factory.mktemp("bank-fixture") / "sequences.sb"
    data = Path(__file__).parent / "testdata"
    bank = SeqBank(path, write=True, engine="lmdb")
    try:
        for name in (
            "NC_024664.1.trunc.fasta",
            "NZ_JAJNFP010000161.1.fasta",
            "NC_036112.1.fasta",
            "NC_044840.1.fasta",
            "NC_010663.1.gb",
            "NC_036113.1.fasta",
        ):
            bank.add_file(data / name)
    finally:
        bank.close()
    return path


@pytest.fixture(autouse=True)
def portable_test_bank(request, monkeypatch):
    if hasattr(request.module, "TEST_BANK_PATH"):
        path = request.getfixturevalue("portable_bank")
        monkeypatch.setattr(request.module, "TEST_BANK_PATH", path)


@pytest.fixture(autouse=True)
def close_stale_banks(monkeypatch):
    """LMDB refuses to open one environment twice in a process; the CLI never closes its banks.

    Close handles left open by earlier code before a bank is reopened, and everything at teardown.
    """
    import seqbank.seqbank as seqbank_module

    opened = []
    real_open_storage = seqbank_module.open_storage

    def open_storage(path, write, requested=None):
        for stale_path, stale in list(opened):
            if stale_path == path.resolve():
                stale.close()
                opened.remove((stale_path, stale))
        storage, engine = real_open_storage(path, write, requested)
        if engine == "lmdb":
            opened.append((path.resolve(), storage))
        return storage, engine

    monkeypatch.setattr(seqbank_module, "open_storage", open_storage)
    yield
    for _, storage in opened:
        try:
            storage.close()
        except Exception:
            pass
