from typing import Union
from functools import cached_property
import numpy as np
import gzip
import time
import tempfile
from pathlib import Path
# import h5py
import zarr
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from zlib import adler32
from attrs import define
from rich.progress import track
import pyfastx
from datetime import datetime

from .transform import seq_to_numpy
from .io import get_file_format, open_path, download_file
from .exceptions import SeqBankError



@define(slots=False)
class SeqBank():
    path:Path
    
    def __getstate__(self):
        # Only returns required elements
        # Needed because h5 files cannot be pickled
        return dict(path=self.path)

    def key(self, accession:str) -> str:
        # Using adler32 for a fast deterministic hash
        accession_hash = str(adler32(accession.encode('ascii')))
        return f"/{accession_hash[-6:-3]}/{accession_hash[-3:]}/{accession}"

    def key_url(self, url:str) -> str:
        accession_hash = str(adler32(url.encode('ascii')))
        filename = Path(url).name
        return f"/urls/{accession_hash[-6:-3]}/{accession_hash[-3:]}/{filename}"

    @cached_property
    def file(self):
        store = zarr.ZipStore(self.path, mode='a')
        return zarr.open(store, mode='a')
        # return h5py.File(self.path, "a", libver='latest')

    def __getitem__(self, accession:str) -> np.ndarray:
        try:
            key = self.key(accession)
            file = self.file
            return file[key][:]
        except Exception as err:
            raise SeqBankError(f"Failed to read {accession} in SeqBank {self.path}:\n{err}")

    def __contains__(self, accession:str) -> bool:
        try:
            return self.key(accession) in self.file
        except Exception:
            return False

    def delete(self, accession:str) -> None:
        key = self.key(accession)
        f = self.file
        if key in f:
            try:
                del f[key]
            except KeyError:
                # f.create_dataset(key, "") # replace data in case there is a problem with it
                # f[key] = "" 
                del f[key]

    def add(self, seq:Union[str, Seq, SeqRecord, np.ndarray], accession:str, all_accessions=None):
        key = self.key(accession)

        if key in self.file:
            return self.file[key]
        
        if isinstance(seq, SeqRecord):
            seq = seq.seq
        if isinstance(seq, Seq):
            seq = str(seq)
        if isinstance(seq, str):
            seq = seq_to_numpy(seq)

        if all_accessions is not None:
            all_accessions.add(accession)
        
        return self.file.create_dataset(
            self.key(accession),
            data=seq,
            dtype="u1",
            # compression="gzip",
            # compression_opts=9,
        )

    def add_file(self, path:Path, format:str=""):
        format = format or get_file_format(path)

        # If fasta or fastq use pyfastx for speed
        if format in ["fasta", "fastq"]:
            for accession, seq in track(pyfastx.Fasta(str(path), build_index=False)):
                self.add(seq, accession)

        with open_path(path) as f:
            for record in track(SeqIO.parse(f, format)):
                self.add(record, record.id)

    def add_url(self, url:str, format:str="", force:bool=False) -> bool:
        url_key = self.key_url(url)
        if url_key in self.file and not force:
            return False
        
        with tempfile.TemporaryDirectory() as tmpdirname:
            local_path = Path(tmpdirname) / Path(url).name
            download_file(url, local_path)
            self.add_file(local_path, format=format)
            
            self.file[url_key] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

            return True

    def get_accessions(self) -> set:
        accessions = set()
        file = self.file
        main_keys = file.keys()
        for key0 in track(main_keys, "Getting accessions"):
            level2_keys = file[f"/{key0}"].keys()
            for key1 in level2_keys:
                dir_accessions = file[f"/{key0}/{key1}"].keys()
                accessions.update(dir_accessions)
        return set(accessions)

    def download_accessions(self, accessions, base_dir:Path, email:str=None, all_accessions=None):
        base_dir.mkdir(exist_ok=True, parents=True)
        local_path = base_dir / "downloaded.fa.gz"

        from Bio import Entrez

        if email:
            Entrez.email = email
        else:
            raise Exception("no email provided")

        print(f"Trying to download '{accessions}'")
        accessions_str = ",".join(accessions)
        try:
            print("trying nucleotide database")
            net_handle = Entrez.efetch(db="nucleotide", id=accessions_str, rettype="fasta", retmode="text")
        except Exception as err:
            print(f'failed {err}')
            print("trying genome database")
            time.sleep(3)
            try:
                net_handle = Entrez.efetch(db="genome", id=accessions_str, rettype="fasta", retmode="text")
            except Exception as err:
                print(f'failed {err}')
                print("trying nuccore database")
                try:
                    net_handle = Entrez.efetch(db="nuccore", id=accessions_str, rettype="fasta", retmode="text")
                except:
                    print(f'failed {err}')
                    return None

        with gzip.open(local_path, "wt") as f:
            f.write(net_handle.read())
        net_handle.close()
        self.add_file(local_path, all_accessions=all_accessions)

    def missing_accessions(self, accessions, get:bool=False):
        missing = set()
        for accession in track(accessions):
            if accession not in self:
                missing.add(accession)
            elif get:
                try:
                    data = self[accession]
                except Exception:
                    missing.add(accession)
        return missing

    def add_accessions(self, accessions, base_dir:Path, email:str=None, all_accessions=None, batch_size=200, sleep:float=1.0):
        to_download = []
        for accession in track(accessions):
            if accession in self:
                continue

            to_download.append(accession)
            if len(to_download) >= batch_size:
                self.download_accessions(to_download, base_dir=base_dir, email=email, all_accessions=all_accessions)
                time.sleep(sleep)
                to_download = []
        
        self.download_accessions(to_download, base_dir=base_dir, email=email, all_accessions=all_accessions)
            # fasta_path = self.individual_accession_path(accession, base_dir=base_dir, email=email)
            # if fasta_path:
            #     self.add_file(fasta_path, all_accessions=all_accessions)

    def individual_accession_path(self, accession: str, base_dir:Path, download: bool = True, email=None) -> Path:
        local_path = base_dir / f"{accession}.fa.gz"
        local_path.parent.mkdir(exist_ok=True, parents=True)
        if download and not local_path.exists():
            from Bio import Entrez

            if email:
                Entrez.email = email
            else:
                raise Exception("no email provided")

            print(f"Trying to download '{accession}'")
            try:
                print("trying nucleotide database")
                net_handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            except Exception as err:
                print(f'failed {err}')
                print("trying genome database")
                time.sleep(1)
                try:
                    net_handle = Entrez.efetch(db="genome", id=accession, rettype="fasta", retmode="text")
                except Exception as err:
                    print(f'failed {err}')
                    print("trying nuccore database")
                    try:
                        net_handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text")
                    except:
                        print(f'failed {err}')
                        return None

            with gzip.open(local_path, "wt") as f:
                f.write(net_handle.read())
            net_handle.close()
        return local_path
