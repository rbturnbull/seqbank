import gzip
from pathlib import Path
import h5py
import tempfile
import shutil

from .io import download_file
from .seqbank import SeqBank

def dfam_url(curated:bool=True, release="current"):
    curated_str = "_curatedonly" if curated else ""
    url = f"https://www.dfam.org/releases/{release}/families/Dfam{curated_str}.h5.gz"
    return url


def add_dfam(seqbank:SeqBank, local_path:Path):
    file = h5py.File(local_path, 'r')
    def visitor_func(name, node):
        if isinstance(node, h5py.Dataset):
            accession = node.attrs['accession']
            seq = node.attrs['consensus']
            seqbank.add(seq=seq, accession=accession)

    file['Families/DF'].visititems(visitor_func)


def download_dfam(seqbank:SeqBank, curated:bool=True, release="current", force:bool=False):
    url = dfam_url(curated, release)
    url_key = seqbank.key_url(url)
    if url_key in seqbank.file and not force:
        print(f"Aleady downloaded: {url}")
        return False
    
    with tempfile.TemporaryDirectory() as tmpdirname:
        local_gzip_path = Path(tmpdirname) / Path(url).name
        try:
            print(f"Downloading Dfam: {url}")
            download_file(url, local_gzip_path)

            # decompress
            local_path = local_gzip_path.with_suffix("")
            with gzip.open(local_gzip_path, 'rb') as infile, open(local_path, 'wb') as outfile:
                shutil.copyfileobj(infile, outfile)

            # Add
            add_dfam(seqbank, local_path)

            seqbank.save_seen_url(url)
        except Exception as err:
            print(f"Failed to add Dfam: {url}: {err}")
            return False

        return True
