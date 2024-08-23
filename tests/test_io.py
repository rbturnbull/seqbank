from seqbank.io import get_file_format, seq_count, download_file
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock
from tempfile import NamedTemporaryFile
from io import BytesIO
import requests

TEST_DATA_PATH = Path(__file__).parent / "testdata"


def test_get_file_format():
    assert get_file_format("test.fasta") == "fasta"
    assert get_file_format("test.fa") == "fasta"
    assert get_file_format("test.fna") == "fasta"    
    assert get_file_format("test.fna.gz") == "fasta"        

    assert get_file_format("test.genbank") == "genbank"        
    assert get_file_format("test.genbank.gz") == "genbank"        
    assert get_file_format("test.gb.gz") == "genbank"        
    assert get_file_format("test.gbk.gz") == "genbank"        

    assert get_file_format("test.embl") == "embl"
    assert get_file_format("test.embl.gz") == "embl"
    assert get_file_format("test.embl.bz2") == "embl"

    assert get_file_format("test.nexus") == "nexus"
    assert get_file_format("test.nxs.gz") == "nexus"

    assert get_file_format("test.tsv.gz") == "tab"        
    
    assert get_file_format("test.fastq") == "fastq"        
    assert get_file_format("test.fq") == "fastq"        

    with pytest.raises(ValueError):
        get_file_format("test.txt")


def test_seq_count():
    seq_count(TEST_DATA_PATH/"NC_010663.1.gb") == 1
    seq_count(TEST_DATA_PATH/"NC_010663.1.gb.bz2") == 1

    seq_count(TEST_DATA_PATH/"5seq.fa.gz") == 5
    seq_count(TEST_DATA_PATH/"5seq.fa") == 5


@pytest.fixture
def dummy_url():
    return "https://example.com/dummyfile"

@pytest.fixture
def dummy_local_path(tmp_path):
    return tmp_path / "dummyfile.gz"

def test_download_file_success(dummy_url, dummy_local_path):
    """Mock the response object from requests.get"""
    mock_response = MagicMock()
    mock_response.iter_content = MagicMock(return_value=[b"dummy data chunk 1", b"dummy data chunk 2"])
    mock_response.__enter__.return_value = mock_response
    mock_response.raise_for_status = MagicMock()
    
    # Patch requests.get to return the mock response
    with patch('seqbank.io.requests.get', return_value=mock_response):
        result_path = download_file(dummy_url, dummy_local_path)
        
        # Check that the file was written correctly
        with open(result_path, 'rb') as f:
            file_content = f.read()
            assert file_content == b"dummy data chunk 1dummy data chunk 2"

        # Check that the function returned the correct path
        assert result_path == dummy_local_path
        
        # Ensure that raise_for_status was called to handle potential HTTP errors
        mock_response.raise_for_status.assert_called_once()

def test_download_file_http_error(dummy_url, dummy_local_path):
    """Mock the response object to raise an HTTPError"""
    mock_response = MagicMock()
    mock_response.__enter__.return_value = mock_response
    mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError("HTTP Error")
    
    with patch('seqbank.io.requests.get', return_value=mock_response):
        with pytest.raises(requests.exceptions.HTTPError, match="HTTP Error"):
            download_file(dummy_url, dummy_local_path)

def test_download_file_no_data(dummy_url, dummy_local_path):
    # Mock the response to return no data
    mock_response = MagicMock()
    mock_response.iter_content = MagicMock(return_value=[])
    mock_response.__enter__.return_value = mock_response
    
    with patch('seqbank.io.requests.get', return_value=mock_response):
        result_path = download_file(dummy_url, dummy_local_path)
        
        # Check that an empty file was created
        with open(result_path, 'rb') as f:
            file_content = f.read()
            assert file_content == b""