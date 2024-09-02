import tempfile
from typer.testing import CliRunner
from seqbank.main import app
from seqbank.seqbank import SeqBank
from unittest.mock import patch, MagicMock
from pathlib import Path
from seqbank.dfam import download_dfam, add_dfam
from plotly.graph_objs import Figure
import shutil
import pytest
import h5py

from .test_seqbank import setup_seqbank

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


def copy_dummy_fasta(url, path):
    if "NC_036113.1.fasta" in url:
        shutil.copy(TEST_DATA_PATH/"NC_036113.1.fasta", path)
    elif "NC_024664.1.trunc.fasta" in url:
        shutil.copy(TEST_DATA_PATH/"NC_024664.1.trunc.fasta", path)


def test_url():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)
        new_path = tmpdirname/"new.sb"

        with patch("seqbank.seqbank.download_file", copy_dummy_fasta):
            result = runner.invoke(app, ["url", str(new_path), "http://example.com/NC_036113.1.fasta"])
            assert result.exit_code == 0

        result = runner.invoke(app, ["ls", str(new_path)])
        assert result.exit_code == 0
        assert result.stdout == '/seqbank/url/http://example.com/NC_036113.1.fasta\nNC_036113.1\n'

def mock_get_refseq_urls():
    return ["http://example.com/NC_036113.1.fasta", "http://example.com/NC_024664.1.trunc.fasta"]

def test_refseq():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)
        new_path = tmpdirname/"new.sb"

        # Patch the get_refseq_urls to return our mock URLs and patch download_file to copy test data
        with patch("seqbank.main.get_refseq_urls", mock_get_refseq_urls):
            with patch("seqbank.seqbank.download_file", copy_dummy_fasta):
                result = runner.invoke(app, ["refseq", str(new_path)])
                assert result.exit_code == 0

        # Check that the sequences have been added correctly
        result = runner.invoke(app, ["ls", str(new_path)])
        assert result.exit_code == 0
        assert 'NC_024664.1\n' in result.stdout
        assert 'NC_036113.1\n' in result.stdout

@pytest.fixture
def seqbank_with_temp_dir():
    """Fixture to create a SeqBank instance with a temporary directory."""
    temp_dir = tempfile.TemporaryDirectory()
    seqbank = SeqBank(path=Path(temp_dir.name), write=True)
    yield seqbank, Path(temp_dir.name)
    temp_dir.cleanup()

@pytest.fixture
def temp_hdf5_file():
    """Fixture to create a temporary HDF5 file with test data."""
    temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.h5')
    file_path = Path(temp_file.name)

    # Create HDF5 file with test data
    with h5py.File(file_path, 'w') as f:
        group = f.create_group('Families/DF')
        dataset = group.create_dataset('test_dataset', data=[], dtype='S1')
        dataset.attrs['accession'] = 'DF000001.1'
        dataset.attrs['consensus'] = 'Mock sequence data for DF000001.1'
        dataset = group.create_dataset('test_dataset2', data=[], dtype='S1')
        dataset.attrs['accession'] = 'DF000002.1'
        dataset.attrs['consensus'] = 'Mock sequence data for DF000002.1'
    
    yield file_path

    # Clean up the temporary file
    file_path.unlink()

def mock_add_dfam(seqbank: SeqBank, local_path: Path):
    """Mock function to simulate adding sequences to SeqBank."""
    seqbank.add = MagicMock()
    seqbank.add(seq="Mock sequence data for DF000001.1", accession="DF000001.1")
    seqbank.add(seq="Mock sequence data for DF000002.1", accession="DF000002.1")

@patch('seqbank.dfam.download_dfam')
def test_dfam(mock_download_dfam, seqbank_with_temp_dir, temp_hdf5_file):
    """Test the dfam command."""
    seqbank, _ = seqbank_with_temp_dir

    # Mock download_dfam to simulate successful download and addition of DFam sequences
    def mock_download_dfam(seqbank, curated=True, release="current", force=False):
        mock_add_dfam(seqbank, temp_hdf5_file)
        return True

    mock_download_dfam.side_effect = mock_download_dfam

    # Run the dfam command
    result = runner.invoke(app, ["dfam", str(seqbank.path)])
    assert result.exit_code == 0, f"DFAM command failed with error: {result.output}"



def test_histogram():
    """Test the save_histogram command."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)

        output_path = tmpdirname/"temp_histogram.png"
        
        runner = CliRunner()

        # Run the save_histogram command
        result = runner.invoke(app, ["histogram", str(TEST_DATA_PATH/"seqbank.sb"), "--output-path", str(output_path)])

        # Check the output message
        assert "Histogram saved to" in result.output

        assert output_path.exists()

        # Assert the command ran successfully
        assert result.exit_code == 0