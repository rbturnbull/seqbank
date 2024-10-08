import typer
from pathlib import Path
from rich.progress import track

from .refseq import get_refseq_urls
from .seqbank import SeqBank
from .dfam import download_dfam

import plotly.graph_objs as go

app = typer.Typer()


@app.command()
def add(path: Path, files: list[Path], format: str = "", filter: Path = None) -> None:
    """Add sequences from a file or list of files to a SeqBank.

    Args:
        path (Path): The path to the SeqBank.
        files (list[Path]): A list of file paths containing sequences.
        format (str, optional): The format of the sequence files. Defaults to "".
        filter (Path, optional): A filter file for sequences. Defaults to None.
    """
    print(f"Opening seqbank '{path}'")
    seqbank = SeqBank(path=path, write=True)
    seqbank.add_files(files, format=format, filter=filter)


@app.command()
def url(path: Path, urls: list[str], format: str = "", max: int = 0, workers: int = -1, tmp_dir: Path = None) -> None:
    """Add sequences from a list of URLs to a SeqBank.

    Args:
        path (Path): The path to the SeqBank.
        urls (list[str]): A list of URLs containing sequences.
        format (str, optional): The format of the sequence files. Defaults to "".
        max (int, optional): Maximum number of sequences to add. Defaults to 0 (all).
        workers (int, optional): Number of workers to use for downloading. Defaults to -1.
        tmp_dir (Path, optional): Temporary directory for downloads. Defaults to None.
    """
    print(f"Opening seqbank '{path}'")
    seqbank = SeqBank(path=path, write=True)

    seqbank.add_urls(urls, format=format, max=max, workers=workers, tmp_dir=tmp_dir)


@app.command()
def delete(path: Path, accessions: list[str]) -> None:
    """Delete sequences from a SeqBank.

    Args:
        path (Path): The path to the SeqBank.
        accessions (list[str]): A list of accessions to delete from the SeqBank.
    """
    print(f"Opening seqbank '{path}'")
    seqbank = SeqBank(path=path, write=True)

    for accession in track(accessions, "Deleting"):
        seqbank.delete(accession)


@app.command()
def refseq(path: Path, max: int = 0, workers: int = -1, tmp_dir: Path = None) -> None:
    """Download all RefSeq sequences to a SeqBank.

    Args:
        path (Path): The path to the SeqBank.
        max (int, optional): Maximum number of sequences to add. Defaults to 0 (all).
        workers (int, optional): Number of workers to use for downloading. Defaults to -1.
        tmp_dir (Path, optional): Temporary directory for downloads. Defaults to None.
    """
    print("Getting RefSeq files list")
    return url(path, get_refseq_urls(tmp_dir=tmp_dir), max=max, workers=workers, tmp_dir=tmp_dir)


@app.command()
def dfam(path: Path, release: str = "current", curated: bool = True) -> bool:
    """Download DFam sequences to a SeqBank.

    Args:
        path (Path): The path to the SeqBank.
        release (str, optional): The DFam release version to download. Defaults to "current".
        curated (bool, optional): Whether to download curated sequences. Defaults to True.

    Returns:
        bool: True if the download and addition were successful, False otherwise.
    """
    print("Getting DFam")
    seqbank = SeqBank(path=path, write=True)
    return download_dfam(seqbank, release=release, curated=curated)


@app.command()
def ls(path: Path) -> None:
    """List accessions in a SeqBank.

    Args:
        path (Path): The path to the SeqBank.
    """
    seqbank = SeqBank(path=path)
    seqbank.ls()


@app.command()
def count(path: Path) -> None:
    """Display the number of accessions in a SeqBank.

    Args:
        path (Path): The path to the SeqBank.
    """
    seqbank = SeqBank(path=path)
    print(len(seqbank))


@app.command()
def cp(path: Path, new: Path) -> None:
    """Copy each sequence from one SeqBank to another.

    Args:
        path (Path): The path to the source SeqBank.
        new (Path): The path to the destination SeqBank.
    """
    print(f"Copying seqbank '{path}' to '{new}'")
    seqbank = SeqBank(path=path)
    new = SeqBank(path=new, write=True)
    seqbank.copy(new)


@app.command()
def export(path: Path, output: Path, format: str = "fasta") -> None:
    """Export a SeqBank to a specified format.

    Args:
        path (Path): The path to the SeqBank.
        output (Path): The path to save the exported sequences.
        format (str, optional): The format for the exported sequences. Defaults to "fasta".
    """
    print(f"Exporting seqbank '{path}' to '{output}' in {format} format")
    seqbank = SeqBank(path=path)
    return seqbank.export(output, format=format)


@app.command()
def histogram(path: Path, output_path: Path = None, show: bool = False, nbins: int = 30) -> None:
    """Generate a histogram of sequence lengths from a SeqBank.

    Args:
        path (Path): The path to the SeqBank.
        output_path (Path, optional): The path to save the histogram. If None, the histogram will be displayed.
        show (bool, optional): Whether to display the histogram. Defaults to False.
        nbins (int, optional): The number of bins for the histogram. Defaults to 30.
    """
    # Load the SeqBank
    seqbank = SeqBank(path=path)

    # Generate the histogram
    fig: go.Figure = seqbank.histogram(nbins=nbins)

    # Save the histogram to the specified output path
    if output_path is None:
        show = True
    else:
        fig.write_image(output_path)
        print(f"Histogram saved to {output_path}")

    if show:
        fig.show()
