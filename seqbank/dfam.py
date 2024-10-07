import gzip
from pathlib import Path
import h5py
import tempfile
import shutil

from .io import download_file
from .seqbank import SeqBank


def dfam_url(curated: bool = True, release: str = "current") -> str:
    """
    Constructs the URL for the Dfam database file.

    Args:
        curated (bool, optional): If True, fetches the curated-only version of the Dfam database. Defaults to True.
        release (str, optional): The Dfam release version. Defaults to "current".

    Returns:
        str: The URL for the Dfam HDF5 file.
    """
    curated_str = "_curatedonly" if curated else ""
    url = f"https://www.dfam.org/releases/{release}/families/Dfam{curated_str}.h5.gz"
    return url


def add_dfam(seqbank: SeqBank, local_path: Path) -> None:
    """
    Adds sequences from a local Dfam HDF5 file to the SeqBank.

    Args:
        seqbank (SeqBank): The SeqBank instance to add sequences to.
        local_path (Path): The path to the local Dfam HDF5 file.
    """
    file = h5py.File(local_path, "r")

    def visitor_func(name, node):
        if isinstance(node, h5py.Dataset):
            accession = node.attrs["accession"]
            seq = node.attrs["consensus"]
            seqbank.add(seq=seq, accession=accession)

    file["Families/DF"].visititems(visitor_func)


def download_dfam(seqbank: SeqBank, curated: bool = True, release: str = "current", force: bool = False) -> bool:
    """
    Downloads the Dfam HDF5 file, decompresses it, and adds sequences to the SeqBank.

    Args:
        seqbank (SeqBank): The SeqBank instance to which the Dfam sequences will be added.
        curated (bool, optional): If True, fetches the curated-only version of the Dfam database. Defaults to True.
        release (str, optional): The Dfam release version. Defaults to "current".
        force (bool, optional): If True, forces the download even if the URL has been previously processed. Defaults to False.

    Returns:
        bool: True if the download and addition were successful, False otherwise.
    """
    url = dfam_url(curated, release)
    url_key = seqbank.key_url(url)

    if url_key in seqbank.file and not force:
        print(f"Already downloaded: {url}")
        return False

    with tempfile.TemporaryDirectory() as tmpdirname:
        local_gzip_path = Path(tmpdirname) / Path(url).name
        try:
            print(f"Downloading Dfam: {url}")
            download_file(url, local_gzip_path)

            # Decompress the gzipped file
            local_path = local_gzip_path.with_suffix("")
            with gzip.open(local_gzip_path, "rb") as infile, open(local_path, "wb") as outfile:
                shutil.copyfileobj(infile, outfile)

            # Add sequences to SeqBank
            add_dfam(seqbank, local_path)

            # Mark URL as processed
            seqbank.save_seen_url(url)
        except Exception as err:
            print(f"Failed to add Dfam: {url}: {err}")
            return False

    return True
