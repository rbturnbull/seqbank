# -*- coding: future_typing -*-

from pathlib import Path
import gzip
import requests

def open_path(path:Path|str):
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def get_file_format(path:Path|str) -> str:    
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix == ".gz":
        suffix = path.suffixes[-2].lower()

    if suffix in [".fa", ".fna", ".fasta"]:
        return "fasta"

    if suffix in [".genbank", ".gb", ".gbk"]:
        return "genbank"

    if suffix in [".embl"]:
        return "embl"

    if suffix in [".tab", ".tsv"]:
        return "tsv"

    if suffix in [".fastq", ".fq"]:
        return "fastq"

    raise ValueError(f"Cannot determine file format of {path}.")


def download_file(url:str, local_path:Path):
    """ Adapted from https://stackoverflow.com/questions/16694907/download-large-file-in-python-with-requests """
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192): 
                f.write(chunk)
    return local_path


