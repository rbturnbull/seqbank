import unittest
import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path
import gzip
import h5py
import tempfile
import shutil

from seqbank.seqbank import SeqBank
from seqbank.dfam import download_dfam, add_dfam, dfam_url
from seqbank.io import download_file

def test_dfam_url_default():
    """Test the default parameters."""
    expected_url = "https://www.dfam.org/releases/current/families/Dfam_curatedonly.h5.gz"
    assert dfam_url() == expected_url


def test_dfam_url_curated_true():
    """Test when curated is True (explicitly set) and release is 'current'."""
    expected_url = "https://www.dfam.org/releases/current/families/Dfam_curatedonly.h5.gz"
    assert dfam_url(curated=True, release="current") == expected_url


def test_dfam_url_curated_false():
    """Test when curated is False."""
    expected_url = "https://www.dfam.org/releases/current/families/Dfam.h5.gz"
    assert dfam_url(curated=False, release="current") == expected_url

def test_dfam_url_custom_release():
    """Test with a custom release value."""
    expected_url = "https://www.dfam.org/releases/2023-07/families/Dfam_curatedonly.h5.gz"
    assert dfam_url(curated=True, release="2023-07") == expected_url

def test_dfam_url_curated_false_custom_release():
    """Test when curated is False with a custom release value."""
    expected_url = "https://www.dfam.org/releases/2023-07/families/Dfam.h5.gz"
    assert dfam_url(curated=False, release="2023-07") == expected_url

@pytest.fixture
def mock_seqbank():
    # Create a mock SeqBank object
    return MagicMock(spec=SeqBank)

@pytest.fixture
def mock_url():
    return "https://www.dfam.org/releases/current/families/Dfam.h5.gz"

@patch('seqbank.dfam.download_file')
@patch('seqbank.dfam.gzip.open')
@patch('seqbank.dfam.h5py.File')
@patch('seqbank.dfam.Path')
def test_download_dfam_already_downloaded(mock_path, mock_h5py, mock_gzip, mock_download_file, mock_seqbank):
    # Setup mocks
    mock_path.return_value = Path("/mock/path/to/Dfam.h5")
    mock_seqbank.key_url.return_value = 'mock_key'
    mock_seqbank.file = {'mock_key': True}  # Simulate file already downloaded
    
    # Call the function
    result = download_dfam(mock_seqbank, curated=True, release="current", force=False)
    
    # Assertions
    assert result is False
    mock_download_file.assert_not_called()
    mock_gzip.assert_not_called()
    mock_h5py.assert_not_called()
    mock_seqbank.add.assert_not_called()
    mock_seqbank.save_seen_url.assert_not_called()

@patch('seqbank.dfam.download_file')
@patch('seqbank.dfam.gzip.open')
@patch('seqbank.dfam.h5py.File')
@patch('seqbank.dfam.Path')
def test_download_dfam_error(mock_path, mock_h5py, mock_gzip, mock_download_file, mock_seqbank):
    # Setup mocks to raise an exception
    mock_download_file.side_effect = Exception("Network error")
    
    # Call the function
    result = download_dfam(mock_seqbank, curated=True, release="current", force=False)
    
    # Assertions
    assert result is False
    mock_download_file.assert_called_once()
    mock_gzip.assert_not_called()
    mock_h5py.assert_not_called()
    mock_seqbank.add.assert_not_called()
    mock_seqbank.save_seen_url.assert_not_called()

@pytest.fixture
def seqbank_with_temp_dir():
    """Fixture to create a SeqBank instance with a temporary directory."""
    temp_dir = tempfile.TemporaryDirectory()
    seqbank = SeqBank(path=Path(temp_dir.name), write=True)
    yield seqbank, Path(temp_dir.name)
    temp_dir.cleanup()

@patch('seqbank.dfam.dfam_url')
@patch('seqbank.dfam.download_file')
@patch('seqbank.dfam.add_dfam')
def test_download_dfam_success(mock_add_dfam, mock_download_file, mock_dfam_url, seqbank_with_temp_dir):
    """Test that the download_dfam function works correctly."""
    seqbank, _ = seqbank_with_temp_dir

    # Mock dfam_url to return a test URL
    mock_dfam_url.return_value = 'http://example.com/test.gz'
    
    # Mock download_file to simulate a successful download
    mock_download_file.return_value = None

    # Mock add_dfam to avoid actual file operations
    mock_add_dfam.return_value = None

    # Mock file creation for the decompressed file
    with patch('gzip.open', unittest.mock.mock_open(read_data=b'Test data')), \
         patch('builtins.open', unittest.mock.mock_open()) as mock_open:
        
        # Run the function
        result = download_dfam(seqbank, curated=True, release="current", force=False)

        # Assertions
        assert result is True
        mock_dfam_url.assert_called_once_with(True, "current")
        mock_download_file.assert_called_once_with('http://example.com/test.gz', unittest.mock.ANY)
        mock_open.assert_called()  # Ensure file operations were attempted
        mock_add_dfam.assert_called()  # Ensure add_dfam was called

@pytest.fixture
def temp_hdf5_file():
    """Fixture to create a temporary HDF5 file with test data."""
    temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.h5')
    file_path = Path(temp_file.name)

    # Create HDF5 file with test data
    with h5py.File(file_path, 'w') as f:
        group = f.create_group('Families/DF')
        dataset = group.create_dataset('test_dataset', data=[], dtype='S1')
        dataset.attrs['accession'] = 'mock_accession'
        dataset.attrs['consensus'] = 'mock_seq'
    
    yield file_path

    # Clean up the temporary file
    file_path.unlink()

def test_add_dfam(temp_hdf5_file):
    """Test the add_dfam function."""
    # Create a mock SeqBank object
    mock_seqbank = MagicMock(spec=SeqBank)
    
    # Call the function
    add_dfam(mock_seqbank, temp_hdf5_file)
    
    # Assertions
    mock_seqbank.add.assert_called_once_with(seq='mock_seq', accession='mock_accession')