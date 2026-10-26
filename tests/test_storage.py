import importlib.util
import json
from concurrent.futures import ThreadPoolExecutor

import pytest
from typer.testing import CliRunner

from seqbank import SeqBank, SeqBankError
from seqbank.main import app
from seqbank.storage import METADATA


@pytest.fixture(params=['rocksdict', 'lmdb', 'speedict'])
def engine(request):
    if importlib.util.find_spec(request.param) is None:
        pytest.skip(f'{request.param} is not installed')
    return request.param


def test_roundtrip(tmp_path, engine):
    path = tmp_path / 'bank.sb'
    bank = SeqBank(path, write=True, engine=engine)
    with ThreadPoolExecutor(4) as pool:
        list(pool.map(lambda i: bank.add('ACGTN', str(i)), range(30)))
    bank.file[bank.key_url('example')] = b'done'
    assert len(bank) == 31
    bank.delete('0')
    bank.close()
    bank = SeqBank(path)
    assert bank['1'] == bytes([1, 2, 3, 4, 0])
    assert '0' not in bank
    assert bank.engine == engine
    assert len(list(bank.items())) == 30
    with pytest.raises(Exception):
        bank.add('A', 'readonly')
    bank.close()
    assert json.loads((path / METADATA).read_text())['engine'] == engine


def test_conversion(tmp_path, engine):
    source = tmp_path / 'source'
    bank = SeqBank(source, write=True, engine=engine)
    bank.add('ACGT', 'sequence')
    bank.file[bank.key_url('example')] = b'done'
    bank.close()
    destination = tmp_path / 'destination'
    target_engine = 'lmdb' if engine != 'lmdb' else 'rocksdict'
    result = CliRunner().invoke(app, ['convert', str(source), str(destination), '--engine', target_engine])
    assert result.exit_code == 0, result.output
    target = SeqBank(destination)
    assert target['sequence'] == bytes([1, 2, 3, 4])
    assert target.file[target.key_url('example')] == b'done'
    target.close()
    result = CliRunner().invoke(app, ['convert', str(source), str(destination), '--engine', target_engine])
    assert result.exit_code != 0


def test_mismatch_and_unknown(tmp_path):
    bank = SeqBank(tmp_path / 'bank', write=True, engine='lmdb')
    bank.add('A', 'a')
    bank.close()
    with pytest.raises(SeqBankError, match='convert'):
        SeqBank(bank.path, engine='rocksdict').file
    with pytest.raises(SeqBankError, match='Unknown engine'):
        SeqBank(tmp_path / 'new', write=True, engine='unknown').file


def test_legacy(tmp_path):
    pytest.importorskip('speedict')
    bank = SeqBank(tmp_path / 'legacy', write=True, engine='speedict')
    bank.add('AC', 'a')
    bank.close()
    (bank.path / METADATA).unlink()
    reopened = SeqBank(bank.path)
    assert reopened['a'] == b'\x01\x02'
    assert reopened.engine == 'speedict'
    reopened.close()
    assert not (bank.path / METADATA).exists()


def test_missing_engine(tmp_path, monkeypatch):
    def missing(name):
        raise ImportError(name)
    monkeypatch.setattr('seqbank.storage.importlib.import_module', missing)
    with pytest.raises(SeqBankError, match="pip install 'seqbank'"):
        SeqBank(tmp_path / 'new', write=True, engine='lmdb').file
    assert not (tmp_path / 'new').exists()


def _process_write(path, number):
    bank = SeqBank(path, write=True)
    for i in range(20):
        bank.add('ACGT', f'{number}-{i}')
    bank.close()


def test_lmdb_multiple_processes(tmp_path):
    import multiprocessing
    path = tmp_path / 'shared'
    bank = SeqBank(path, write=True, engine='lmdb')
    bank.file
    bank.close()
    context = multiprocessing.get_context('spawn')
    processes = [context.Process(target=_process_write, args=(path, i)) for i in range(3)]
    for process in processes:
        process.start()
    for process in processes:
        process.join(30)
        if process.is_alive():
            process.terminate()
            process.join()
            pytest.fail('Writer process timed out')
        assert process.exitcode == 0
    reader = SeqBank(path)
    assert len(reader) == 60
    reader.close()


def test_bad_metadata(tmp_path):
    path = tmp_path / 'bank'
    path.mkdir()
    (path / METADATA).write_text('{"engine": "lmdb", "format_version": 99}')
    with pytest.raises(SeqBankError, match='Invalid SeqBank metadata'):
        SeqBank(path).file


def test_parallel_import(tmp_path, engine):
    files = []
    for i in range(4):
        path = tmp_path / f'{i}.fasta'
        path.write_text(''.join(f'>{i}-{j}\nACGTN\n' for j in range(1100)))
        files.append(path)
    bank = SeqBank(tmp_path / 'imported', write=True, engine=engine)
    bank.add_files(files, workers=4)
    assert len(bank) == 4400
    assert bank['3-1099'] == b'\x01\x02\x03\x04\x00'
    bank.close()


def test_existing_legacy_fixture():
    from pathlib import Path
    pytest.importorskip('speedict')
    bank = SeqBank(Path(__file__).parent / 'testdata' / 'seqbank.sb')
    try:
        assert len(bank) == 6
        assert bank.engine == 'speedict'
        assert len(bank['NC_024664.1']) == 840
    finally:
        bank.close()


def test_engine_enum(tmp_path):
    from seqbank import Engine
    bank = SeqBank(tmp_path / 'enum-bank', write=True, engine=Engine.LMDB)
    bank.add('AC', 'a')
    bank.close()
    reader = SeqBank(bank.path)
    assert reader['a'] == b'\x01\x02'
    assert reader.engine is Engine.LMDB
    reader.close()
    assert json.loads((bank.path / METADATA).read_text())['engine'] == 'lmdb'


@pytest.mark.parametrize('command', ['add', 'add-sequence-from-file', 'url', 'refseq', 'dfam', 'convert'])
def test_cli_engine_choices(command):
    result = CliRunner().invoke(app, [command, '--help'], color=False, terminal_width=160, env={'COLUMNS': '160'})
    assert result.exit_code == 0
    for name in ('speedict', 'rocksdict', 'lmdb'):
        assert name in result.output


def test_cli_rejects_invalid_engine_before_opening(tmp_path):
    source = tmp_path / 'missing'
    destination = tmp_path / 'destination'
    result = CliRunner().invoke(app, ['convert', str(source), str(destination), '--engine', 'invalid'])
    assert result.exit_code == 2
    assert 'Invalid value' in result.output
    assert not destination.exists()


def test_long_url_tracking_and_conversion(tmp_path, engine):
    url = 'https://example.com/' + 'x' * 700
    other_url = url + 'different'
    bank = SeqBank(tmp_path / 'long-urls', write=True, engine=engine)
    assert not bank.seen_url(url)
    bank.save_seen_url(url)
    assert bank.seen_url(url)
    assert not bank.seen_url(other_url)
    assert all(len(key) <= 511 for key in bank.file.keys())
    bank.close()
    bank = SeqBank(bank.path)
    assert bank.seen_url(url)
    target = SeqBank(tmp_path / 'converted-urls', write=True,
                     engine='rocksdict' if engine == 'lmdb' else 'lmdb')
    bank.copy(target)
    assert target.seen_url(url)
    assert not target.seen_url(other_url)
    del target.file[target.key_url(url)]
    assert not target.seen_url(url)
    bank.close()
    target.close()


def test_legacy_raw_long_url_conversion(tmp_path):
    url = 'https://example.com/' + 'x' * 700
    source = SeqBank(tmp_path / 'old-rocks', write=True, engine='rocksdict')
    # Simulate an unmodified long tracking key from a pre-release bank.
    source.file.db[source.key_url(url)] = b'old'
    assert source.seen_url(url)
    source.file[source.key_url(url)] = b'updated'
    assert source.file[source.key_url(url)] == b'updated'
    target = SeqBank(tmp_path / 'new-lmdb', write=True, engine='lmdb')
    source.copy(target)
    assert target.seen_url(url)
    assert target.file[target.key_url(url)] == b'updated'
    source.close()
    target.close()
