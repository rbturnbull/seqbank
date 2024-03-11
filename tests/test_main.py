import tempfile
from typer.testing import CliRunner
from seqbank.main import app
from unittest.mock import patch
from pathlib import Path
import shutil

TEST_DATA_PATH = Path(__file__).parent / "testdata"

runner = CliRunner()


def test_export():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)

        exported_path = tmpdirname/"test_main_export.fasta"
        result = runner.invoke(app, ["export", str(TEST_DATA_PATH/"seqbank.sb"), str(exported_path)])
        assert result.exit_code == 0
        assert exported_path.exists()
        text = exported_path.read_text()
        assert text.count(">") == 6
        assert ">NC_024664.1\nTAAAAAGAAAAAAATTTT" in text
        assert ">NZ_JAJNFP010000161.1\nTAATATTTGCTTTCATTTCTAAATAG" in text
        assert ">NC_010663.1\nAGTTTTAAACCTCTGATCGAAC" in text
        assert ">NC_036112.1\nGAAATACCCAATATC" in text
        assert ">NC_036113.1\nGAAATACCCAATATCTTGTTCC" in text
        assert ">NC_044840.1\nGGGCGAACGACGGGAATTGAACCCGC" in text
        assert len(text) > 5000


def test_count():
    result = runner.invoke(app, ["count", str(TEST_DATA_PATH/"seqbank.sb")])
    assert result.exit_code == 0
    assert "6\n" in result.stdout


def test_ls():
    result = runner.invoke(app, ["ls", str(TEST_DATA_PATH/"seqbank.sb")])
    assert result.exit_code == 0
    assert result.stdout == 'NC_010663.1\nNC_024664.1\nNC_036112.1\nNC_036113.1\nNC_044840.1\nNZ_JAJNFP010000161.1\n'


def test_cp():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)

        new_path = tmpdirname/"new.sb"
        assert new_path.exists() == False
        result = runner.invoke(app, ["cp", str(TEST_DATA_PATH/"seqbank.sb"), str(new_path)])
        assert result.exit_code == 0
        assert new_path.exists()

        result = runner.invoke(app, ["ls", str(new_path)])
        assert result.exit_code == 0
        assert result.stdout == 'NC_010663.1\nNC_024664.1\nNC_036112.1\nNC_036113.1\nNC_044840.1\nNZ_JAJNFP010000161.1\n'


def test_delete():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)
        new_path = tmpdirname/"new.sb"

        shutil.copytree(TEST_DATA_PATH/"seqbank.sb", new_path)

        result = runner.invoke(app, ["delete", str(new_path), "NC_010663.1"])
        assert result.exit_code == 0

        result = runner.invoke(app, ["ls", str(new_path)])
        assert result.exit_code == 0
        assert "NC_010663.1\n" not in result.stdout
        assert 'NC_024664.1\nNC_036112.1\nNC_036113.1\nNC_044840.1\nNZ_JAJNFP010000161.1\n' in result.stdout


def test_add():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)
        new_path = tmpdirname/"new.sb"

        result = runner.invoke(app, ["add", str(new_path), str(TEST_DATA_PATH/"NC_036113.1.fasta"), str(TEST_DATA_PATH/"NC_024664.1.trunc.fasta")])
        assert result.exit_code == 0

        result = runner.invoke(app, ["count", str(new_path)])
        assert result.exit_code == 0
        assert "2\n" in result.stdout


def test_add_fail():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)
        new_path = tmpdirname/"new.sb"

        result = runner.invoke(app, ["add", str(new_path), str(TEST_DATA_PATH/"NC_036113.1.NOT_HERE"), str(TEST_DATA_PATH/"NC_024664.1.trunc.fasta")])
        assert result.exit_code == 1


def copy_NC_036113(url, path):
    shutil.copy(TEST_DATA_PATH/"NC_036113.1.fasta", path)


def test_url():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)
        new_path = tmpdirname/"new.sb"

        with patch("seqbank.seqbank.download_file", copy_NC_036113):
            result = runner.invoke(app, ["url", str(new_path), "http://example.com/NC_036113.1.fasta"])
            assert result.exit_code == 0

        result = runner.invoke(app, ["ls", str(new_path)])
        assert result.exit_code == 0
        assert result.stdout == '/seqbank/url/http://example.com/NC_036113.1.fasta\nNC_036113.1\n'


        
