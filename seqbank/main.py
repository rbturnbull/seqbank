
import typer
from pathlib import Path
from typing import List
from rich.progress import track

from .refseq import get_refseq_urls
from .seqbank import SeqBank
from .dfam import download_dfam

app = typer.Typer()

@app.command()
def add(path:Path, files:List[Path], format:str="", filter:Path=None):
    """ Add sequences from a file or list of files to a seqbank """
    print(f"Opening seqbank '{path}'")
    seqbank = SeqBank(path=path, write=True)
    seqbank.add_files(files, format=format, filter=filter)


@app.command()
def url(path:Path, urls:List[str], format:str="", max:int=0, workers:int=-1):
    """ Add sequences from a URL to a seqbank """
    print(f"Opening seqbank '{path}'")
    seqbank = SeqBank(path=path, write=True)
    
    seqbank.add_urls(urls, format=format, max=max, workers=workers)


@app.command()
def delete(path:Path, accessions:List[str]):
    """ Deletes sequences from a seqbank """

    print(f"Opening seqbank '{path}'")
    seqbank = SeqBank(path=path, write=True)
    
    for accession in track(accessions, "Deleting"):
        seqbank.delete(accession)


@app.command()
def refseq(path:Path, max:int=0, workers:int=-1):
    """ Download all RefSeq sequences to a seqbank """    
    print("Getting RefSeq files list")
    return url(path, get_refseq_urls(), max=max, workers=workers)


@app.command()
def dfam(path:Path, release:str="current", curated:bool=True):
    """ Download all RefSeq sequences to a seqbank """    
    print("Getting DFam")
    seqbank = SeqBank(path=path, write=True)
    return download_dfam(seqbank, release=release, curated=curated)


@app.command()
def ls(path:Path):
    """ List accessions in a seqbank """  
    seqbank = SeqBank(path=path)
    seqbank.ls()


@app.command()
def count(path:Path):
    """ Displays the number of accessions in a seqbank """  
    seqbank = SeqBank(path=path)
    print(len(seqbank))


@app.command()
def cp(path:Path, new:Path):
    """ Copies each sequence from one seqbank to another. """
    print(f"Copying seqbank '{path}' to '{new}'")
    seqbank = SeqBank(path=path)
    new = SeqBank(path=new, write=True)
    seqbank.copy(new)


@app.command()
def export(path:Path, output:Path, format:str="fasta"):
    print(f"Exporting seqbank '{path}' to '{output}' in {format} format")
    seqbank = SeqBank(path=path)
    return seqbank.export(output, format=format)
