# -*- coding: future_typing -*-

from typing import Union, List
from functools import cached_property
import numpy as np
import gzip
import time
import tempfile
from pathlib import Path
from joblib import Parallel, delayed
# import zarr
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from attrs import define
from rich.progress import track
import pyfastx
from rich.progress import Progress, TimeElapsedColumn, MofNCompleteColumn
from datetime import datetime
from .transform import seq_to_bytes, bytes_to_str
from .io import get_file_format, open_path, download_file
from .exceptions import SeqBankError
# from rocksdict import Rdict, Options
from speedict import Rdict, Options, DBCompressionType, AccessType
import atexit


@define(slots=False)
class SeqBank():
    path:Path
    write:bool = False

    def __attrs_post_init__(self):
        if not self.write and not self.path.exists():
            raise FileNotFoundError(f"Cannot find SeqBank file at path: {self.path}")
    
    def __getstate__(self):
        # Only returns required elements
        # Needed because h5 files cannot be pickled
        return dict(path=self.path)

    def key(self, accession:str) -> str:
        return bytes(accession, "ascii")

    def key_url(self, url:str) -> str:
        return self.key("/seqbank/url/"+url)

    def close(self):
        print(f"Closing SeqBank '{self.path}'")
        try:
            self._db.close()
        except Exception:
            pass

    @cached_property
    def file(self):
        options = Options(raw_mode=True)
        options.set_compression_type(DBCompressionType.none())
        # options.set_cache_index_and_filter_blocks(True)
        options.set_optimize_filters_for_hits(True)
        options.optimize_for_point_lookup(1024)
        options.set_max_open_files(500)

        self._db = Rdict(
            path=str(self.path), 
            options=options, 
            access_type=AccessType.read_write() if self.write else AccessType.read_only()
        )

        atexit.register(self.close)
        return self._db
        # store = zarr.DBMStore(self.path, open=dbm.gnu.open)
        # store = zarr.ZipStore(self.path, mode='a')
        # return zarr.open(store, mode='a')
        # return h5py.File(self.path, "a", libver='latest')

    def __len__(self):
        count = 0
        for _ in track(self.file.keys()):
            count += 1
        
        return count

    def __getitem__(self, accession:str) -> np.ndarray:
        try:
            key = self.key(accession)
            file = self.file
            return file[key]
            # return np.frombuffer(file[key], dtype="u1")
        except Exception as err:
            raise SeqBankError(f"Failed to read {accession} in SeqBank {self.path}:\n{err}")

    def items(self):
        for k,v in self.file.items():
            yield k, np.frombuffer(v, dtype="u1")

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

    def add(self, seq:Union[str, Seq, SeqRecord, np.ndarray], accession:str, all_accessions=None) -> None:
        key = self.key(accession)
        # if key in self.file:
        #     return
        
        if isinstance(seq, SeqRecord):
            seq = seq.seq
        if isinstance(seq, Seq):
            seq = str(seq)
        if isinstance(seq, str):
            seq = seq_to_bytes(seq)
        
        # seq = compress(seq, self.compression)
        self.file[key] = seq

    def add_file(self, path:Path, format:str="", progress=None, overall_task=None) -> None:
        format = format or get_file_format(path)
        progress = progress or Progress()

        # If fasta or fastq use pyfastx for speed
        if format in ["fasta", "fastq"]:
            total = sum(1 for _ in pyfastx.Fasta(str(path), build_index=False))
            task = progress.add_task(f"[magenta]{path.name}", total=total)

            for accession, seq in pyfastx.Fasta(str(path), build_index=False):
                self.add(seq, accession)
                progress.update(task, advance=1)
        else:
            with open_path(path) as f:
                total = sum(1 for _ in SeqIO.parse(f, format))
            
            with open_path(path) as f:
                task = progress.add_task(f"[magenta]{path.name}", total=total)
                for record in SeqIO.parse(f, format):
                    self.add(record, record.id)
                    progress.update(task, advance=1)

        progress.update(task, visible=False)
        if overall_task is not None:
            progress.update(overall_task, advance=1)

        print('Added', path.name)

    def seen_url(self, url:str) -> bool:
        return self.key_url(url) in self.file

    def save_seen_url(self, url:str):
        url_key = self.key_url(url)
        self.file[url_key] = bytes(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "ascii")

    def add_url(self, url:str, progress=None, format:str="", force:bool=False, overall_task=None) -> bool:
        url_key = self.key_url(url)
        if url_key in self.file and not force:
            return False
        
        with tempfile.TemporaryDirectory() as tmpdirname:
            local_path = Path(tmpdirname) / Path(url).name
            try:
                download_file(url, local_path)
                self.add_file(local_path, format=format, progress=progress, overall_task=overall_task)
                self.save_seen_url(url)
            except Exception as err:
                print(f"Failed to add URL: {url}: {err}")
                return False

            return True

    def get_accessions(self) -> set:
        accessions = set()
        file = self.file
        main_keys = file.keys()
        for key in file.keys():
            accession = key.decode("ascii")
            if not accession.startswith("/seqbank/"):
                accessions.update([accession])
        return accessions

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

    def missing(self, accessions, get:bool=False) -> set:
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

    def add_urls(self, urls:List[str], max:int=0, format:str="", force:bool=False, workers:int=-1):
        # only add the URLs that haven't been seen before
        urls_to_add = []
        for url in urls:
            if not self.seen_url(url):
                urls_to_add.append(url)

            # truncate URLs list to `max` if requested
            if max and len(urls_to_add) >= max:
                break

        with Progress(*Progress.get_default_columns(), TimeElapsedColumn(), MofNCompleteColumn()) as progress:
            parallel = Parallel(n_jobs=workers, prefer="threads")
            add_url = delayed(self.add_url)
            overall_task = progress.add_task(f"[bold red]Adding URLs", total=len(urls_to_add))
            parallel(add_url(url, progress=progress, format=format, force=force, overall_task=overall_task) for url in urls_to_add)

    def ls(self):
        breakpoint()

    def add_files(self, files:List[str], max:int=0, format:str="", workers:int=1):
        with Progress(*Progress.get_default_columns(), TimeElapsedColumn(), MofNCompleteColumn()) as progress:
            parallel = Parallel(n_jobs=workers, prefer="threads")
            add_file = delayed(self.add_file)
            overall_task = progress.add_task(f"[bold red]Adding files", total=len(files))
            parallel(add_file(file, progress=progress, format=format, overall_task=overall_task) for file in files)

    def copy(self, other):
        for k,v in track(self.file.items()):
            other.file[k] = v

    def numpy(self, accession) -> np.ndarray:
        """
        Returns the seq data for an accession as an unsigned char NumPy array.
        """
        return np.frombuffer(self[accession], dtype="u1")

    def string(self, accession:str) -> str:
        """
        Returns the seq data for an accession as a string.
        """
        data = self[accession]
        return bytes_to_str(data)

    def record(self, accession:str) -> SeqRecord:
        """
        Returns a BioPython SeqRecord object for an accession.
        """
        record = SeqRecord(
            Seq(self.string(accession)),
            id=accession,
            description="",
        )
        return record

    def export(self, output:Path|str, format:str="", accessions:List[str]|str|Path|None=None):
        """
        Exports the data in a seqbank to a file using BioPython
        """
        accessions = accessions or self.get_accessions()

        # Read list of accessions if given file
        if isinstance(accessions, (str,Path)):
            accessions = Path(accessions).read_text().strip().split("\n")

        format = format or get_file_format(output)
        with open(output, "w") as f:
            for accession in accessions:
                SeqIO.write(self.record(accession), f, format)

