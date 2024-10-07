from pathlib import Path
import gzip
import bz2
from Bio import SeqIO
import requests
import tempfile


def open_path(path: Path | str):
    """
    Opens a file, optionally decompressing it based on its extension.

    Args:
        path (Path | str): The path to the file, which can be a compressed (.gz, .bz2) or uncompressed file.

    Returns:
        file object: A file object opened for reading.

    Raises:
        ValueError: If the file has an unsupported compression extension.
    """
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix == ".gz":
        return gzip.open(path, "rt")

    if suffix == ".bz2":
        return bz2.open(path, "rt")

    return open(path, "rt")


def get_file_format(path: Path | str) -> str:
    """
    Determines the sequence file format based on the file extension.

    Supported file formats and their extensions:
        - FASTA: .fa, .fna, .fasta
        - GenBank: .genbank, .gb, .gbk
        - EMBL: .embl
        - Tab-delimited: .tab, .tsv
        - Nexus: .nexus, .nxs
        - FASTQ: .fastq, .fq

    Args:
        path (Path | str): The path to the file.

    Returns:
        str: The detected file format (e.g., 'fasta', 'genbank', 'embl', etc.).

    Raises:
        ValueError: If the file format cannot be determined.
    """
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


def seq_count(path: Path | str) -> int:
    """
    Counts the number of sequences in a file.

    Args:
        path (Path | str): The path to the file.

    Returns:
        int: The total number of sequences in the file.
    """
    file_format = get_file_format(path)
    with open_path(path) as f:
        total = sum(1 for _ in SeqIO.parse(f, file_format))
    return total


def download_file(url: str, local_path: Path) -> Path:
    """
    Downloads a file from a URL and saves it to a specified local path.

    Args:
        url (str): The URL of the file to download.
        local_path (Path): The path where the file should be saved.

    Returns:
        Path: The path to the downloaded file.

    Raises:
        requests.exceptions.HTTPError: If the download fails.
    """
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return local_path


class TemporaryDirectory(tempfile.TemporaryDirectory):
    """
    A custom TemporaryDirectory class that can optionally create a parent directory if a prefix is provided.

    Args:
        prefix (str | Path | None): Optional prefix for the temporary directory.
    """

    def __init__(self, prefix: str | Path | None = None, *args, **kwargs):
        self._created_dirs = []  # Track directories that were created
        if isinstance(prefix, Path):
            # Resolve path to string and ensure parent directories are created
            prefix.mkdir(parents=True, exist_ok=True)
            self._created_dirs.append(prefix)  # Track created parent directories

            # Ensure prefix ends with a "/"
            prefix = str(prefix.resolve().absolute())
            if not prefix.endswith("/"):
                prefix += "/"

        super().__init__(prefix=prefix, *args, **kwargs)

    def cleanup(self):
        """
        Cleans up the temporary directory and removes any created parent directories if they are empty.
        """
        super().cleanup()
        for created_dir in self._created_dirs:
            if not any(created_dir.iterdir()):  # Only remove if the directory is empty
                created_dir.rmdir()

    def __enter__(self) -> Path:
        """
        Enters the context of the temporary directory.

        Returns:
            Path: The path to the temporary directory.
        """
        return Path(super().__enter__())
