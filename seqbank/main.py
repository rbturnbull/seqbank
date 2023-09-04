
import typer
from pathlib import Path
from typing import List
from rich.progress import track

from .refseq import get_refseq_urls
from .seqbank import SeqBank

app = typer.Typer()

@app.command()
def add(path:Path, files:List[Path], format:str=""):
    """ Add sequences from a filepath to a seqbank """
    print(f"Opening seqbank '{path}'")
    seqbank = SeqBank(path=path)
    
    for path in files:
        print("Adding file", path)
        seqbank.add_file(path, format=format)


@app.command()
def url(path:Path, urls:List[str], format:str="", max:int=0):
    """ Add sequences from a URL to a seqbank """
    print(f"Opening seqbank '{path}'")
    seqbank = SeqBank(path=path)
    
    count = 0
    for url in urls:
        print("Adding URL", url)
        result = seqbank.add_url(url, format=format)
        count += int(result)

        if count >= max:
            print(f"Maximum number of URLs reached: {max}")
            break


@app.command()
def delete(path:Path, accessions:List[str]):
    """ Deletes sequences from a seqbank """

    print(f"Opening seqbank '{path}'")
    seqbank = SeqBank(path=path)
    
    for accesion in track(accessions, "Deleting"):
        seqbank.delete(accesion)


@app.command()
def refseq(path:Path, max:int=0):
    """ Download all RefSeq sequences to a seqbank """    
    print("Getting RefSeq files list")
    return url(path, get_refseq_urls(), max=max)



