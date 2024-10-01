from pathlib import Path
import re
from .io import download_file, TemporaryDirectory

def get_refseq_filenames(tmp_dir: str | Path | None = None) -> list[str]:
    """
    Retrieves a list of RefSeq genomic filenames from the NCBI FTP site.

    Args:
        tmp_dir (str | Path | None, optional): 
            The directory to create a temporary folder in. If None, a default temporary directory is used.
    
    Returns:
        list[str]: A list of filenames sorted numerically by version.
    """
    with TemporaryDirectory(prefix=tmp_dir) as dirname:
        local_path = dirname / 'refseq_complete.html'

        download_file("https://ftp.ncbi.nlm.nih.gov/refseq/release/complete/", local_path)
        text = local_path.read_text()

        # Find all .genomic.fna.gz filenames
        filenames = re.findall(r'>(.*?.genomic.fna.gz)<\/a>', text)
        
        # Sort filenames numerically by version number (second part of filename)
        filenames = sorted(filenames, key=lambda filename: int(filename.split(".")[1]))

        return filenames
    
    
def get_refseq_urls(tmp_dir: str | Path | None = None) -> list[str]:
    """
    Retrieves a list of URLs for RefSeq genomic files from the NCBI FTP site.

    Args:
        tmp_dir (str | Path | None, optional): 
            The directory to create a temporary folder in. If None, a default temporary directory is used.

    Returns:
        list[str]: A list of URLs for the RefSeq genomic files.
    """
    filenames = get_refseq_filenames(tmp_dir=tmp_dir)
    return [f"https://ftp.ncbi.nlm.nih.gov/refseq/release/complete/{filename}" for filename in filenames]