import typer
from typing import Annotated
from pathlib import Path
from rich.progress import track
import plotly.io as pio

from .refseq import get_refseq_urls
from .seqbank import SeqBank
from .storage import Engine
from .dfam import download_dfam

pio.defaults.mathjax = None


app = typer.Typer(pretty_exceptions_enable=False)


@app.command()
def add(
    path: Annotated[
        Path, typer.Argument(help="SeqBank database directory; existing banks detect their engine automatically.")
    ],
    files: Annotated[
        list[Path], typer.Argument(help="Sequence files to import (FASTA, FASTQ, or another supported format).")
    ],
    format: Annotated[str, typer.Option(help="Input sequence format; infer from the filename when omitted.")] = "",
    filter: Annotated[
        Path, typer.Option(help="File listing accessions to include; omit to import all sequences.")
    ] = None,
    engine: Annotated[
        Engine,
        typer.Option(
            help="Storage engine for a new bank (default: rocksdict if installed, otherwise lmdb); must match an existing bank."
        ),
    ] = None,
) -> None:
    """Add sequences from a file or list of files to a SeqBank."""
    print(f"Opening seqbank '{path}'")
    seqbank = SeqBank(path=path, write=True, engine=engine)
    seqbank.add_files(files, format=format, filter=filter)


@app.command()
def add_sequence_from_file(
    path: Annotated[
        Path, typer.Argument(help="SeqBank database directory; existing banks detect their engine automatically.")
    ],
    accession: Annotated[str, typer.Argument(help="Accession to assign to the imported sequence.")],
    file: Annotated[Path, typer.Argument(help="Input file containing exactly one sequence.")],
    format: Annotated[str, typer.Option(help="Input sequence format; infer from the filename when omitted.")] = "",
    engine: Annotated[
        Engine,
        typer.Option(
            help="Storage engine for a new bank (default: rocksdict if installed, otherwise lmdb); must match an existing bank."
        ),
    ] = None,
) -> None:
    """Add a sequence from a file with a single sequence to a SeqBank."""
    print(f"Opening seqbank '{path}'")
    seqbank = SeqBank(path=path, write=True, engine=engine)
    seqbank.add_sequence_from_file(accession, file, format=format)
    print(f"Added {accession} from {file.name}")


@app.command()
def url(
    path: Annotated[
        Path, typer.Argument(help="SeqBank database directory; existing banks detect their engine automatically.")
    ],
    urls: Annotated[list[str], typer.Argument(help="URLs of sequence files to download and import.")],
    format: Annotated[str, typer.Option(help="Input sequence format; infer from the filename when omitted.")] = "",
    max: Annotated[int, typer.Option(help="Maximum number of URLs to import; 0 imports all URLs.")] = 0,
    workers: Annotated[
        int, typer.Option(help="Number of download/import threads; -1 uses all available CPU cores.")
    ] = -1,
    tmp_dir: Annotated[
        Path, typer.Option(help="Directory for temporary downloads; omit to use the system temporary directory.")
    ] = None,
    engine: Annotated[
        Engine,
        typer.Option(
            help="Storage engine for a new bank (default: rocksdict if installed, otherwise lmdb); must match an existing bank."
        ),
    ] = None,
) -> None:
    """Add sequences from a list of URLs to a SeqBank."""
    print(f"Opening seqbank '{path}'")
    seqbank = SeqBank(path=path, write=True, engine=engine)

    seqbank.add_urls(urls, format=format, max=max, workers=workers, tmp_dir=tmp_dir)


@app.command()
def delete(
    path: Annotated[
        Path, typer.Argument(help="SeqBank database directory; existing banks detect their engine automatically.")
    ],
    accessions: Annotated[list[str], typer.Argument(help="Accessions to delete from the bank.")],
) -> None:
    """Delete sequences from a SeqBank."""
    print(f"Opening seqbank '{path}'")
    seqbank = SeqBank(path=path, write=True)

    for accession in track(accessions, "Deleting"):
        seqbank.delete(accession)


@app.command()
def refseq(
    path: Annotated[
        Path, typer.Argument(help="SeqBank database directory; existing banks detect their engine automatically.")
    ],
    max: Annotated[int, typer.Option(help="Maximum number of URLs to import; 0 imports all URLs.")] = 0,
    workers: Annotated[
        int, typer.Option(help="Number of download/import threads; -1 uses all available CPU cores.")
    ] = -1,
    tmp_dir: Annotated[
        Path, typer.Option(help="Directory for temporary downloads; omit to use the system temporary directory.")
    ] = None,
    engine: Annotated[
        Engine,
        typer.Option(
            help="Storage engine for a new bank (default: rocksdict if installed, otherwise lmdb); must match an existing bank."
        ),
    ] = None,
) -> None:
    """Download all RefSeq sequences to a SeqBank."""
    print("Getting RefSeq files list")
    return url(path, get_refseq_urls(tmp_dir=tmp_dir), max=max, workers=workers, tmp_dir=tmp_dir, engine=engine)


@app.command()
def dfam(
    path: Annotated[
        Path, typer.Argument(help="SeqBank database directory; existing banks detect their engine automatically.")
    ],
    release: Annotated[str, typer.Option(help="DFam release to download, such as current.")] = "current",
    curated: Annotated[
        bool, typer.Option(help="Download only curated DFam families; disable to include all families.")
    ] = True,
    engine: Annotated[
        Engine,
        typer.Option(
            help="Storage engine for a new bank (default: rocksdict if installed, otherwise lmdb); must match an existing bank."
        ),
    ] = None,
) -> bool:
    """Download DFam sequences to a SeqBank."""
    print("Getting DFam")
    seqbank = SeqBank(path=path, write=True, engine=engine)
    return download_dfam(seqbank, release=release, curated=curated)


@app.command()
def ls(
    path: Annotated[
        Path, typer.Argument(help="SeqBank database directory; existing banks detect their engine automatically.")
    ],
) -> None:
    """List accessions in a SeqBank."""
    seqbank = SeqBank(path=path)
    seqbank.ls()


@app.command()
def count(
    path: Annotated[
        Path, typer.Argument(help="SeqBank database directory; existing banks detect their engine automatically.")
    ],
) -> None:
    """Display the number of accessions in a SeqBank."""
    seqbank = SeqBank(path=path)
    print(len(seqbank))


@app.command()
def merge(
    path: Annotated[
        Path, typer.Argument(help="SeqBank database directory; existing banks detect their engine automatically.")
    ],
    new: Annotated[
        Path,
        typer.Argument(
            help="Destination bank; new banks use RocksDict if installed (otherwise LMDB), existing banks retain their engine."
        ),
    ],
) -> None:
    """Merge every sequence from one SeqBank into another."""
    print(f"Merging seqbank '{path}' into '{new}'")
    seqbank = SeqBank(path=path)
    new = SeqBank(path=new, write=True)
    seqbank.copy(new)


@app.command()
def export(
    path: Annotated[
        Path, typer.Argument(help="SeqBank database directory; existing banks detect their engine automatically.")
    ],
    output: Annotated[Path, typer.Argument(help="Destination file for exported sequences.")],
    format: Annotated[str, typer.Option(help="Output sequence format, such as fasta or fastq.")] = "fasta",
) -> None:
    """Export a SeqBank to a specified format."""
    print(f"Exporting seqbank '{path}' to '{output}' in {format} format")
    seqbank = SeqBank(path=path)
    return seqbank.export(output, format=format)


@app.command()
def histogram(
    path: Annotated[
        Path, typer.Argument(help="SeqBank database directory; existing banks detect their engine automatically.")
    ],
    output_path: Annotated[
        Path, typer.Option(help="Save the histogram image here (requires Chrome); omit to display it.")
    ] = None,
    show: Annotated[bool, typer.Option(help="Display the histogram, including when an image is also saved.")] = False,
    nbins: Annotated[int, typer.Option(help="Requested number of histogram bins.")] = 30,
    min: Annotated[
        int, typer.Option(help="Minimum sequence length to include, in bases; 0 disables the lower bound.")
    ] = 0,
    max: Annotated[
        int, typer.Option(help="Maximum sequence length to include, in bases; 0 disables the upper bound.")
    ] = 0,
) -> None:
    """Generate a histogram of sequence lengths from a SeqBank."""
    # Load the SeqBank
    seqbank = SeqBank(path=path)

    # Generate the histogram
    fig = seqbank.histogram(nbins=nbins, min=min, max=max)

    # Save the histogram to the specified output path
    if output_path is None:
        show = True
    else:
        fig.write_image(output_path)
        print(f"Histogram saved to {output_path}")

    if show:
        fig.show()


@app.command()
def convert(
    path: Annotated[
        Path, typer.Argument(help="SeqBank database directory; existing banks detect their engine automatically.")
    ],
    destination: Annotated[Path, typer.Argument(help="New destination database directory; must not already exist.")],
    engine: Annotated[
        Engine,
        typer.Option(
            help="Storage engine for the destination; RocksDict and Speedict require the rocksdict and speedict extras."
        ),
    ],
) -> None:
    """Convert all records into a new directory using the selected engine."""
    if destination.exists():
        raise typer.BadParameter("Destination must not exist")
    source = SeqBank(path)
    target = SeqBank(destination, write=True, engine=engine)
    try:
        source.copy(target)
    finally:
        source.close()
        target.close()
