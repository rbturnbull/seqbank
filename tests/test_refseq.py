from unittest.mock import patch
import pytest
from pathlib import Path

from seqbank.refseq import get_refseq_filenames, get_refseq_urls


def mock_download_refseq_listing(url:str, path:Path):
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
    assert url == "https://ftp.ncbi.nlm.nih.gov/refseq/release/complete/"
    path.write_text(html_content)


@patch('seqbank.refseq.download_file', mock_download_refseq_listing)
def test_get_refseq_filenames():    
    # Test the get_refseq_filenames function
    filenames = get_refseq_filenames(tmp_dir=None)
    assert filenames == ['refseq.2.genomic.fna.gz', 'refseq.10.genomic.fna.gz', 'refseq.100.genomic.fna.gz']


@patch('seqbank.refseq.download_file', mock_download_refseq_listing)
def test_get_refseq_urls():
    urls = get_refseq_urls(tmp_dir=None)
    expected_urls = [
        'https://ftp.ncbi.nlm.nih.gov/refseq/release/complete/refseq.2.genomic.fna.gz',
        'https://ftp.ncbi.nlm.nih.gov/refseq/release/complete/refseq.10.genomic.fna.gz',
        'https://ftp.ncbi.nlm.nih.gov/refseq/release/complete/refseq.100.genomic.fna.gz',
    ]
    assert urls == expected_urls