from unittest.mock import patch
import pytest
from pathlib import Path

from seqbank.refseq import get_refseq_filenames, get_refseq_urls

@pytest.fixture
def mock_download_file():
    with patch('seqbank.io.download_file') as mock_download:
        yield mock_download

@patch('pathlib.Path.read_text')
def test_get_refseq_filenames(mock_read_text, mock_download_file):
    # Mocking the HTML content of the refseq_complete.html file
    html_content = """
    <html>
    <body>
    <a href="refseq.10.genomic.fna.gz">refseq.10.genomic.fna.gz</a>
    <a href="refseq.2.genomic.fna.gz">refseq.2.genomic.fna.gz</a>
    <a href="refseq.100.genomic.fna.gz">refseq.100.genomic.fna.gz</a>
    </body>
    </html>
    """
    
    # Mock the read_text method to return the HTML content
    mock_read_text.return_value = html_content
    
    # Test the get_refseq_filenames function
    filenames = get_refseq_filenames(tmp_dir=None)
    assert filenames == ['refseq.2.genomic.fna.gz', 'refseq.10.genomic.fna.gz', 'refseq.100.genomic.fna.gz']

def test_get_refseq_urls(mock_download_file):
    # Mocking the output of get_refseq_filenames
    with patch('seqbank.refseq.get_refseq_filenames', return_value=['refseq.2.genomic.fna.gz', 'refseq.10.genomic.fna.gz']):
        urls = get_refseq_urls(tmp_dir=None)
        expected_urls = [
            'https://ftp.ncbi.nlm.nih.gov/refseq/release/complete/refseq.2.genomic.fna.gz',
            'https://ftp.ncbi.nlm.nih.gov/refseq/release/complete/refseq.10.genomic.fna.gz'
        ]
        assert urls == expected_urls