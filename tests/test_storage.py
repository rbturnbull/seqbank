import importlib.util
import json
from concurrent.futures import ThreadPoolExecutor

import pytest
from typer.testing import CliRunner

from seqbank import SeqBank, SeqBankError
from seqbank.main import app
from seqbank.storage import METADATA


@pytest.fixture(params=["rocksdict", "lmdb", "speedict"])
def engine(request):
    if importlib.util.find_spec(request.param) is None:
        pytest.skip(f"{request.param} is not installed")
    return request.param


def test_roundtrip(tmp_path, engine):
    path = tmp_path / "bank.sb"
    bank = SeqBank(path, write=True, engine=engine)
    with ThreadPoolExecutor(4) as pool:
        list(pool.map(lambda i: bank.add("ACGTN", str(i)), range(30)))
    bank.file[bank.key_url("example")] = b"done"
    assert len(bank) == 31
    bank.delete("0")
    bank.close()
    bank = SeqBank(path)
    assert bank["1"] == bytes([1, 2, 3, 4, 0])
    assert "0" not in bank
    assert bank.engine == engine
    assert len(list(bank.items())) == 30
    with pytest.raises(Exception):
        bank.add("A", "readonly")
    bank.close()
    assert json.loads((path / METADATA).read_text())["engine"] == engine


def test_conversion(tmp_path, engine):
    source = tmp_path / "source"
    bank = SeqBank(source, write=True, engine=engine)
    bank.add("ACGT", "sequence")
    bank.file[bank.key_url("example")] = b"done"
    bank.close()
    destination = tmp_path / "destination"
    if engine == "lmdb":
        pytest.importorskip("rocksdict")
    target_engine = "lmdb" if engine != "lmdb" else "rocksdict"
    result = CliRunner().invoke(app, ["convert", str(source), str(destination), "--engine", target_engine])
    assert result.exit_code == 0, result.output
    target = SeqBank(destination)
    assert target["sequence"] == bytes([1, 2, 3, 4])
    assert target.file[target.key_url("example")] == b"done"
    target.close()
    result = CliRunner().invoke(app, ["convert", str(source), str(destination), "--engine", target_engine])
    assert result.exit_code != 0


def test_mismatch_and_unknown(tmp_path):
    bank = SeqBank(tmp_path / "bank", write=True, engine="lmdb")
    bank.add("A", "a")
    bank.close()
    with pytest.raises(SeqBankError, match="convert"):
        SeqBank(bank.path, engine="rocksdict").file
    with pytest.raises(SeqBankError, match="Unknown engine"):
        SeqBank(tmp_path / "new", write=True, engine="unknown").file


def test_legacy(tmp_path):
    pytest.importorskip("speedict")
    bank = SeqBank(tmp_path / "legacy", write=True, engine="speedict")
    bank.add("AC", "a")
    bank.close()
    (bank.path / METADATA).unlink()
    reopened = SeqBank(bank.path)
    assert reopened["a"] == b"\x01\x02"
    assert reopened.engine == "speedict"
    reopened.close()
    assert not (bank.path / METADATA).exists()


def test_missing_engine(tmp_path, monkeypatch):
    def missing(name):
        raise ImportError(name)

    monkeypatch.setattr("seqbank.storage.importlib.import_module", missing)
    with pytest.raises(SeqBankError, match="pip install seqbank"):
        SeqBank(tmp_path / "new", write=True, engine="lmdb").file
    assert not (tmp_path / "new").exists()


def _process_write(path, number):
    bank = SeqBank(path, write=True)
    for i in range(20):
        bank.add("ACGT", f"{number}-{i}")
    bank.close()


def test_lmdb_multiple_processes(tmp_path):
    import multiprocessing

    path = tmp_path / "shared"
    bank = SeqBank(path, write=True, engine="lmdb")
    bank.file
    bank.close()
    context = multiprocessing.get_context("spawn")
    processes = [context.Process(target=_process_write, args=(path, i)) for i in range(3)]
    for process in processes:
        process.start()
    for process in processes:
        process.join(30)
        if process.is_alive():
            process.terminate()
            process.join()
            pytest.fail("Writer process timed out")
        assert process.exitcode == 0
    reader = SeqBank(path)
    assert len(reader) == 60
    reader.close()


def test_bad_metadata(tmp_path):
    path = tmp_path / "bank"
    path.mkdir()
    (path / METADATA).write_text('{"engine": "lmdb", "format_version": 99}')
    with pytest.raises(SeqBankError, match="Invalid SeqBank metadata"):
        SeqBank(path).file


def test_parallel_import(tmp_path, engine):
    files = []
    for i in range(4):
        path = tmp_path / f"{i}.fasta"
        path.write_text("".join(f">{i}-{j}\nACGTN\n" for j in range(1100)))
        files.append(path)
    bank = SeqBank(tmp_path / "imported", write=True, engine=engine)
    bank.add_files(files, workers=4)
    assert len(bank) == 4400
    assert bank["3-1099"] == b"\x01\x02\x03\x04\x00"
    bank.close()


def test_existing_legacy_fixture():
    from pathlib import Path

    pytest.importorskip("speedict")
    bank = SeqBank(Path(__file__).parent / "testdata" / "seqbank.sb")
    try:
        assert len(bank) == 6
        assert bank.engine == "speedict"
        assert len(bank["NC_024664.1"]) == 840
    finally:
        bank.close()


def test_engine_enum(tmp_path):
    from seqbank import Engine

    bank = SeqBank(tmp_path / "enum-bank", write=True, engine=Engine.LMDB)
    bank.add("AC", "a")
    bank.close()
    reader = SeqBank(bank.path)
    assert reader["a"] == b"\x01\x02"
    assert reader.engine is Engine.LMDB
    reader.close()
    assert json.loads((bank.path / METADATA).read_text())["engine"] == "lmdb"


@pytest.mark.parametrize("command", ["add", "add-sequence-from-file", "url", "refseq", "dfam", "convert"])
def test_cli_engine_choices(command):
    result = CliRunner().invoke(app, [command, "--help"], color=False, terminal_width=160, env={"COLUMNS": "160"})
    assert result.exit_code == 0
    for name in ("speedict", "rocksdict", "lmdb"):
        assert name in result.output


def test_cli_rejects_invalid_engine_before_opening(tmp_path):
    source = tmp_path / "missing"
    destination = tmp_path / "destination"
    result = CliRunner().invoke(app, ["convert", str(source), str(destination), "--engine", "invalid"])
    assert result.exit_code == 2
    assert "Invalid value" in result.output
    assert not destination.exists()


def test_long_url_tracking_and_conversion(tmp_path, engine):
    url = "https://example.com/" + "x" * 700
    other_url = url + "different"
    bank = SeqBank(tmp_path / "long-urls", write=True, engine=engine)
    assert not bank.seen_url(url)
    bank.save_seen_url(url)
    assert bank.seen_url(url)
    assert not bank.seen_url(other_url)
    assert all(len(key) <= 511 for key in bank.file.keys())
    bank.close()
    bank = SeqBank(bank.path)
    assert bank.seen_url(url)
    if engine == "lmdb":
        pytest.importorskip("rocksdict")
    target = SeqBank(tmp_path / "converted-urls", write=True, engine="rocksdict" if engine == "lmdb" else "lmdb")
    bank.copy(target)
    assert target.seen_url(url)
    assert not target.seen_url(other_url)
    del target.file[target.key_url(url)]
    assert not target.seen_url(url)
    bank.close()
    target.close()


def test_legacy_raw_long_url_conversion(tmp_path):
    pytest.importorskip("rocksdict")
    url = "https://example.com/" + "x" * 700
    source = SeqBank(tmp_path / "old-rocks", write=True, engine="rocksdict")
    # Simulate an unmodified long tracking key from a pre-release bank.
    source.file.db[source.key_url(url)] = b"old"
    assert source.seen_url(url)
    source.file[source.key_url(url)] = b"updated"
    assert source.file[source.key_url(url)] == b"updated"
    target = SeqBank(tmp_path / "new-lmdb", write=True, engine="lmdb")
    source.copy(target)
    assert target.seen_url(url)
    assert target.file[target.key_url(url)] == b"updated"
    source.close()
    target.close()


@pytest.fixture
def no_rocksdict(monkeypatch):
    """Make ``import rocksdict`` fail, as in an environment without the extra."""
    real_import = importlib.import_module

    def fake_import(name, *args, **kwargs):
        if name == "rocksdict":
            raise ImportError(name)
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr("seqbank.storage.importlib.import_module", fake_import)


def test_default_engine(monkeypatch):
    from seqbank.storage import Engine, default_engine

    if importlib.util.find_spec("rocksdict") is not None:
        assert default_engine() is Engine.ROCKSDICT
    monkeypatch.setitem(__import__("sys").modules, "rocksdict", None)
    assert default_engine() is Engine.LMDB


def test_new_bank_falls_back_to_lmdb(tmp_path, no_rocksdict):
    path = tmp_path / "bank"
    bank = SeqBank(path, write=True)
    bank.add("ACGTN", "a")
    assert bank.engine == "lmdb"
    bank.close()
    assert json.loads((path / METADATA).read_text())["engine"] == "lmdb"
    reader = SeqBank(path)
    assert reader["a"] == bytes([1, 2, 3, 4, 0])
    reader.close()


def test_existing_rocksdict_bank_without_rocksdict(tmp_path, no_rocksdict):
    path = tmp_path / "bank"
    path.mkdir()
    (path / METADATA).write_text('{"engine": "rocksdict", "format_version": 1}')
    with pytest.raises(SeqBankError, match=r"seqbank\[rocksdict\]"):
        SeqBank(path).file


def test_requested_rocksdict_without_rocksdict(tmp_path, no_rocksdict):
    with pytest.raises(SeqBankError, match=r"seqbank\[rocksdict\]"):
        SeqBank(tmp_path / "new", write=True, engine="rocksdict").file
    assert not (tmp_path / "new").exists()


def test_speedict_hint(tmp_path, monkeypatch):
    def missing(name):
        raise ImportError(name)

    monkeypatch.setattr("seqbank.storage.importlib.import_module", missing)
    with pytest.raises(SeqBankError, match=r"seqbank\[speedict\].*Python <3.13"):
        SeqBank(tmp_path / "new", write=True, engine="speedict").file


def test_cli_lmdb_only(tmp_path, no_rocksdict):
    runner = CliRunner()
    bank = tmp_path / "bank.sb"
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">one\nACGT\n>two\nTTGA\n")
    result = runner.invoke(app, ["add", str(bank), str(fasta)])
    assert result.exit_code == 0, result.output
    assert json.loads((bank / METADATA).read_text())["engine"] == "lmdb"
    result = runner.invoke(app, ["ls", str(bank)])
    assert result.exit_code == 0 and result.stdout == "one\ntwo\n"
    result = runner.invoke(app, ["count", str(bank)])
    assert result.exit_code == 0 and "2" in result.stdout
    out = tmp_path / "out.fasta"
    result = runner.invoke(app, ["export", str(bank), str(out)])
    assert result.exit_code == 0, result.output
    assert ">one\nACGT" in out.read_text()


def test_unidentifiable_existing_directory(tmp_path):
    path = tmp_path / "bank"
    path.mkdir()
    (path / "something").write_text("not a bank")
    with pytest.raises(SeqBankError, match="Cannot identify"):
        SeqBank(path).file


def test_legacy_detected_without_metadata(tmp_path):
    # A directory with a RocksDB CURRENT file but no metadata is a pre-metadata Speedict bank.
    path = tmp_path / "bank"
    path.mkdir()
    (path / "CURRENT").write_text("MANIFEST-000001\n")
    with pytest.raises(SeqBankError, match="Bank uses speedict, not lmdb"):
        SeqBank(path, engine="lmdb").file


def test_storage_missing_keys(tmp_path, engine):
    bank = SeqBank(tmp_path / "bank", write=True, engine=engine)
    with pytest.raises(KeyError):
        bank.file[b"missing"]
    bank.close()


def test_lmdb_delete_missing_key(tmp_path):
    bank = SeqBank(tmp_path / "bank", write=True, engine="lmdb")
    with pytest.raises(KeyError):
        del bank.file[b"missing"]
    bank.close()


def test_long_key_lookup_by_raw_key(tmp_path, engine):
    from seqbank.storage import _storage_key

    bank = SeqBank(tmp_path / "bank", write=True, engine=engine)
    key = bank.key_url("https://example.com/" + "x" * 700)
    bank.file[key] = b"value"
    assert _storage_key(key) != key
    assert bank.file[key] == b"value"
    bank.close()


def test_metadata_created_concurrently(tmp_path, monkeypatch):
    """Another writer creating the metadata file first is not an error."""
    from pathlib import Path

    real_open = Path.open

    def racing_open(self, mode="r", *args, **kwargs):
        if self.name == METADATA and mode == "x":
            # Another process wins the race between the existence check and creation.
            self.write_text('{"engine": "lmdb", "format_version": 1}')
        return real_open(self, mode, *args, **kwargs)

    monkeypatch.setattr(Path, "open", racing_open)
    bank = SeqBank(tmp_path / "bank", write=True, engine="lmdb")
    bank.add("A", "a")
    bank.close()
    assert json.loads((bank.path / METADATA).read_text())["engine"] == "lmdb"


def test_metadata_write_failure_closes_storage(tmp_path, monkeypatch):
    import seqbank.storage as storage

    closed = []
    real_storage = storage.LMDBStorage

    class TrackingStorage(real_storage):
        def close(self):
            closed.append(True)
            super().close()

    def failing_dump(*args, **kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(storage, "LMDBStorage", TrackingStorage)
    monkeypatch.setattr(storage.json, "dump", failing_dump)
    with pytest.raises(OSError, match="disk full"):
        storage.open_storage(tmp_path / "bank", write=True, requested="lmdb")
    assert closed == [True]


def test_copy_in_batches(tmp_path, engine):
    source = SeqBank(tmp_path / "source", write=True, engine=engine)
    source.file.put_many((f"/seqbank/{i}".encode(), b"\x01\x02") for i in range(2500))
    destination = SeqBank(tmp_path / "destination", write=True, engine="lmdb")
    source.copy(destination)
    assert len(destination) == 2500
    source.close()
    destination.close()


def test_copy_empty_bank_creates_destination(tmp_path):
    source = SeqBank(tmp_path / "source", write=True, engine="lmdb")
    source.file
    destination_path = tmp_path / "destination"
    destination = SeqBank(destination_path, write=True, engine="lmdb")
    source.copy(destination)
    assert (destination_path / METADATA).exists()
    assert len(destination) == 0
    source.close()
    destination.close()
