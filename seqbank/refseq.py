import tempfile
from pathlib import Path
import re

from .io import download_file

def get_refseq_filenames():
    with tempfile.TemporaryDirectory() as dirname:
        local_path = Path(dirname) / 'refseq_complete.html'

        download_file("https://ftp.ncbi.nlm.nih.gov/refseq/release/complete/", local_path)
        text = local_path.read_text()

        filenames = re.findall(r'>(.*?.genomic.fna.gz)<\/a>', text)
        filenames = sorted(filenames, key=lambda filename:int(filename.split(".")[1]))

        return filenames
    
    
def get_refseq_urls():
    return [f"https://ftp.ncbi.nlm.nih.gov/refseq/release/complete/{filename}" for filename in get_refseq_filenames()]
