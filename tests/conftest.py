"""Build a portable fixture; retain separate tests for legacy Speedict banks."""
from pathlib import Path

import pytest

from seqbank import SeqBank


@pytest.fixture(scope='session')
def portable_bank(tmp_path_factory):
    path = tmp_path_factory.mktemp('bank-fixture') / 'sequences.sb'
    data = Path(__file__).parent / 'testdata'
    bank = SeqBank(path, write=True, engine='rocksdict')
    try:
        for name in (
            'NC_024664.1.trunc.fasta', 'NZ_JAJNFP010000161.1.fasta',
            'NC_036112.1.fasta', 'NC_044840.1.fasta', 'NC_010663.1.gb',
            'NC_036113.1.fasta',
        ):
            bank.add_file(data / name)
    finally:
        bank.close()
    return path


@pytest.fixture(autouse=True)
def portable_test_bank(request, monkeypatch):
    if hasattr(request.module, 'TEST_BANK_PATH'):
        path = request.getfixturevalue('portable_bank')
        monkeypatch.setattr(request.module, 'TEST_BANK_PATH', path)
