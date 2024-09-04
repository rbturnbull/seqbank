
from pathlib import Path
import gzip
import bz2
from Bio import SeqIO
import requests
import tempfile

def open_path(path:Path|str):
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix == ".gz":
        return gzip.open(path, "rt")
    
    if suffix == ".bz2":
        return bz2.open(path, "rt")
    
    return open(path, "rt")


def get_file_format(path:Path|str) -> str:    
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix in [".gz", ".bz2"]:
        suffix = path.suffixes[-2].lower()

    if suffix in [".fa", ".fna", ".fasta"]:
        return "fasta"

    if suffix in [".genbank", ".gb", ".gbk"]:
        return "genbank"

    if suffix in [".embl"]:
        return "embl"

    if suffix in [".tab", ".tsv"]:
        return "tab"

    if suffix in [".nexus", ".nxs"]:
        return "nexus"

    if suffix in [".fastq", ".fq"]:
        return "fastq"

    raise ValueError(f"Cannot determine file format of {path}.")


def seq_count(path:Path|str) -> int:
    format = get_file_format(path)
    with open_path(path) as f:
        total = sum(1 for _ in SeqIO.parse(f, format))
    return total


def download_file(url:str, local_path:Path):
    """ Adapted from https://stackoverflow.com/questions/16694907/download-large-file-in-python-with-requests """
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192): 
                f.write(chunk)
    return local_path


class TemporaryDirectory(tempfile.TemporaryDirectory):
    def __init__(self, prefix:str|Path|None=None, *args, **kwargs):
        if isinstance(prefix, Path):
            # resolve path to string
            prefix = str(prefix.resolve().absolute())
            
        super().__init__(prefix=prefix, *args, **kwargs)
    