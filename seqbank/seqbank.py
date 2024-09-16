
from typing import Union, List, Set
from functools import cached_property
import numpy as np
import gzip
import time
from pathlib import Path
from joblib import Parallel, delayed
import plotly.express as px
import plotly.graph_objs as go
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from attrs import define
from rich.progress import track
import pyfastx
from rich.progress import Progress, TimeElapsedColumn, MofNCompleteColumn
from datetime import datetime
from speedict import Rdict, Options, DBCompressionType, AccessType
import atexit

from .transform import seq_to_bytes, bytes_to_str
from .io import get_file_format, open_path, download_file, seq_count, TemporaryDirectory
from .exceptions import SeqBankError
from .utils import parse_filter



@define(slots=False)
class SeqBank():
    path:Path
    write:bool = False

    def __attrs_post_init__(self):
        self.path = Path(self.path).expanduser()
        if not self.write and not self.path.exists():
            raise FileNotFoundError(f"Cannot find SeqBank file at path: {self.path}")

    def key(self, accession:str) -> str:
        return bytes(accession, "ascii")

    def key_url(self, url:str) -> str:
        return self.key("/seqbank/url/"+url)

    def close(self):
        try:
            self._db.close()
        except Exception:
            pass

    @cached_property
    def file(self):
        options = Options(raw_mode=True)
        options.set_compression_type(DBCompressionType.none())
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
            del f[key]

    def add(self, seq:Union[str, Seq, SeqRecord, np.ndarray], accession:str) -> None:
        key = self.key(accession)
        
        if isinstance(seq, SeqRecord):
            seq = seq.seq
        if isinstance(seq, Seq):
            seq = str(seq)
        if isinstance(seq, str):
            seq = seq_to_bytes(seq)
        
        self.file[key] = seq

    def add_file(self, path:Path, format:str="", progress=None, overall_task=None, filter:Path|list|set|None=None) -> None:
        filter = parse_filter(filter)
        format = format or get_file_format(path)
        progress = progress or Progress()

        # If fasta or fastq use pyfastx for speed
        if format in ["fasta", "fastq"]:
            total = sum(1 for _ in pyfastx.Fasta(str(path), build_index=False))
            task = progress.add_task(f"[magenta]{path.name}", total=total)

            for accession, seq in pyfastx.Fasta(str(path), build_index=False):
                if filter and accession not in filter:
                    continue
                self.add(seq, accession)
                progress.update(task, advance=1)
        else:
            total = seq_count(path)
            
            with open_path(path) as f:
                task = progress.add_task(f"[magenta]{path.name}", total=total)
                for record in SeqIO.parse(f, format):
                    if filter and record.id not in filter:
                        continue

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

    def add_url(self, url:str, progress=None, format:str="", force:bool=False, overall_task=None, tmp_dir:str|Path|None=None) -> bool:
        url_key = self.key_url(url)
        if url_key in self.file and not force:
            return False
        
        with TemporaryDirectory(prefix=tmp_dir) as tmpdirname:
            local_path = tmpdirname / Path(url).name
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
        for key in file.keys():
            accession = key.decode("ascii")
            if not accession.startswith("/seqbank/"):
                accessions.update([accession])
        return accessions

    def missing(self, accessions) -> set:
        missing = set()
        for accession in track(accessions):
            if accession not in self:
                missing.add(accession)
        return missing

    def add_urls(self, urls:List[str], max:int=0, format:str="", force:bool=False, workers:int=-1, tmp_dir:str|Path|None=None):
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
            parallel(add_url(url, progress=progress, format=format, force=force, overall_task=overall_task, tmp_dir=tmp_dir) for url in urls_to_add)

    def ls(self):
        for k in self.file.keys():
            print(k.decode("ascii"))

    def add_files(self, files:List[str], max:int=0, format:str="", workers:int=1, filter:Path|list|set|None=None):
        filter = parse_filter(filter)

        with Progress(*Progress.get_default_columns(), TimeElapsedColumn(), MofNCompleteColumn()) as progress:
            parallel = Parallel(n_jobs=workers, prefer="threads")
            add_file = delayed(self.add_file)
            overall_task = progress.add_task(f"[bold red]Adding files", total=len(files))
            parallel(add_file(file, progress=progress, format=format, overall_task=overall_task, filter=filter) for file in files)

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

    def lengths_dict(self) -> dict[str, int]:
        """
        Returns a dictionary where the keys are the accessions and the values 
        are the corresponding lengths of each sequence.
        """
        accession_lengths = {}

        for accession in self.get_accessions():
            # Retrieve the sequence
            sequence = self[accession]
            # Store the length of the sequence in the dictionary
            accession_lengths[accession] = len(sequence)

        return accession_lengths
    
    def histogram(self, nbins: int = 30) -> go.Figure:
        """
        Creates a histogram of the lengths of all sequences and returns the Plotly figure object.
        """
        # Get the dictionary of accession lengths
        accession_lengths = self.lengths_dict()

        # Extract the lengths from the dictionary
        lengths = list(accession_lengths.values())

        # Create the histogram using Plotly Express
        fig = px.histogram(lengths, nbins=nbins, title="Histogram of Sequence Lengths")
        
        # Add labels and customize the layout, removing the legend
        fig.update_layout(
            xaxis_title="Sequence Length",
            yaxis_title="Count",
            showlegend=False  # Remove the legend
        )

        return fig