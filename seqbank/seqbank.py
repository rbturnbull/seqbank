
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

    def __attrs_post_init__(self) -> None:
        """Initializes the SeqBank object after attributes are set.

        Expands the user path for the SeqBank file and checks if it exists when not in write mode.

        Raises:
            FileNotFoundError: If the SeqBank file is not found at the given path when write mode is disabled.
        """
        self.path = Path(self.path).expanduser()
        if not self.write and not self.path.exists():
            raise FileNotFoundError(f"Cannot find SeqBank file at path: {self.path}")

    def key(self, accession: str) -> str:
        """Generates a key for a given accession.

        Args:
            accession (str): Accession string used to generate the key.

        Returns:
            str: A byte-encoded key as a string.
        """
        return bytes(accession, "ascii")


    def key_url(self, url: str) -> str:
        """Generates a key for a given URL.

        Args:
            url (str): The URL for which the key is generated.

        Returns:
            str: A byte-encoded key string prefixed with '/seqbank/url/'.
        """
        return self.key("/seqbank/url/" + url)


    def close(self) -> None:
        """Closes the SeqBank database connection.

        Attempts to close the database connection, if it's open, and silently handles any exceptions.
        """
        try:
            self._db.close()
        except Exception:
            pass

    @cached_property
    def file(self) -> Rdict:
        """Initializes and configures the Rdict database for sequence storage.

        Configures options for the database, such as compression type, optimization, and maximum open files.
        Registers the close method to be executed upon program exit to ensure the database is closed.

        Returns:
            Rdict: The configured Rdict database object.
        """
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

    def __len__(self) -> int:
        """Calculates the total number of items in the SeqBank.

        Iterates over all the keys in the database and counts them.

        Returns:
            int: The number of entries in the SeqBank.
        """
        count = 0
        for _ in track(self.file.keys()):
            count += 1
        return count

    def __getitem__(self, accession: str) -> np.ndarray:
        """Retrieves the sequence data associated with a given accession.

        Args:
            accession (str): The accession key to look up in the SeqBank.

        Returns:
            np.ndarray: The sequence data stored in the SeqBank for the given accession.

        Raises:
            SeqBankError: If the accession cannot be read or an error occurs during retrieval.
        """
        try:
            key = self.key(accession)
            file = self.file
            return file[key]
        except Exception as err:
            raise SeqBankError(f"Failed to read {accession} in SeqBank {self.path}:\n{err}")

    def items(self):
        """Yields all key-value pairs from the SeqBank.

        The key is the accession, and the value is the sequence data in NumPy array format.

        Yields:
            tuple: A tuple containing the accession (key) and the sequence data (value) as a NumPy array.
        """
        for k, v in self.file.items():
            yield k, np.frombuffer(v, dtype="u1")

    def __contains__(self, accession: str) -> bool:
        """Checks if a given accession exists in the SeqBank.

        Args:
            accession (str): The accession to check for existence.

        Returns:
            bool: True if the accession exists in the SeqBank, otherwise False.
        """
        try:
            return self.key(accession) in self.file
        except Exception:
            return False

    def delete(self, accession: str) -> None:
        """Deletes a sequence entry from the SeqBank by its accession.

        Args:
            accession (str): The accession of the sequence to delete.

        Returns:
            None
        """
        key = self.key(accession)
        f = self.file
        if key in f:
            del f[key]

    def add(self, seq: Union[str, Seq, SeqRecord, np.ndarray], accession: str) -> None:
        """Adds a sequence to the SeqBank with a given accession.

        This method accepts a sequence in various formats (string, Seq, SeqRecord, or NumPy array)
        and stores it in the SeqBank after appropriate conversion to byte format.

        Args:
            seq (Union[str, Seq, SeqRecord, np.ndarray]): The sequence to add to the SeqBank. 
                It can be a string, Bio.Seq object, SeqRecord, or a NumPy array.
            accession (str): The accession key for the sequence to be stored under.

        Returns:
            None
        """
        key = self.key(accession)
        
        if isinstance(seq, SeqRecord):
            seq = seq.seq
        if isinstance(seq, Seq):
            seq = str(seq)
        if isinstance(seq, str):
            seq = seq_to_bytes(seq)
        
        self.file[key] = seq

    def add_file(self, path:Path, format:str="", progress=None, overall_task=None, filter:Path|list|set|None=None) -> None:
        """Adds sequences from a file to the SeqBank.

        This method processes a sequence file in various formats (e.g., FASTA, FASTQ), optionally filtering specific accessions, 
        and adds the sequences to the SeqBank. Progress tracking is available for large file imports.

        Args:
            path (Path): The path to the sequence file.
            format (str, optional): The format of the sequence file (e.g., "fasta", "fastq"). If not provided, it will be auto-detected.
            progress (Progress, optional): A rich progress bar to display the import progress. Defaults to None.
            overall_task (Union[int, None], optional): An optional task ID for tracking the overall progress. Defaults to None.
            filter (Union[Path, list, set, None], optional): A filter for selecting specific accessions. Defaults to None.

        Returns:
            None
        """
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

    def seen_url(self, url: str) -> bool:
        """Checks if a given URL has been seen (i.e., processed) before and present in the SeqBank.

        Args:
            url (str): The URL to check.

        Returns:
            bool: True if the URL has been seen (i.e., exists in the SeqBank), otherwise False.
        """
        return self.key_url(url) in self.file

    def save_seen_url(self, url: str) -> None:
        """Saves a URL as 'seen' by adding it to the SeqBank with a timestamp.

        Args:
            url (str): The URL to save as seen.

        Returns:
            None
        """
        url_key = self.key_url(url)
        self.file[url_key] = bytes(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "ascii")

    def add_url(self, url:str, progress=None, format:str="", force:bool=False, overall_task=None, tmp_dir:str|Path|None=None) -> bool:
        """Downloads and adds sequences from a URL to the SeqBank.

        This method downloads a file from a given URL, processes it to extract sequences, and adds them to the SeqBank.
        If the URL has already been processed, it can be skipped unless `force=True` is provided.

        Args:
            url (str): The URL to download the sequence file from.
            progress (Progress, optional): A rich progress bar to display progress of the download and file processing. Defaults to None.
            format (str, optional): The format of the sequence file (e.g., "fasta", "fastq"). If not provided, it will be auto-detected. Defaults to "".
            force (bool, optional): Whether to force downloading and processing the URL even if it has been seen before. Defaults to False.
            overall_task (Union[int, None], optional): An optional task ID for tracking overall progress. Defaults to None.
            tmp_dir (Union[str, Path, None], optional): A temporary directory to store the downloaded file. Defaults to None.

        Returns:
            bool: True if the URL was successfully processed and added, False otherwise.
        """
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

    def get_accessions(self) -> Set[str]:
        """Retrieves all accessions stored in the SeqBank.

        This method iterates through the SeqBank database keys and collects all accessions that do not belong 
        to the internal '/seqbank/' namespace.

        Returns:
            Set[str]: A set of all accessions present in the SeqBank.
        """
        accessions = set()
        file = self.file

        for key in file.keys():
            accession = key.decode("ascii")
            if not accession.startswith("/seqbank/"):
                accessions.update([accession])

        return accessions

    def missing(self, accessions: Union[List[str], Set[str]]) -> Set[str]:
        """Finds accessions that are not present in the SeqBank.

        This method checks a list or set of accessions and returns those that are missing from the SeqBank.

        Args:
            accessions (Union[List[str], Set[str]]): A list or set of accessions to check for presence in the SeqBank.

        Returns:
            Set[str]: A set of accessions that are missing from the SeqBank.
        """
        missing = set()

        for accession in track(accessions):
            if accession not in self:
                missing.add(accession)

        return missing

    def add_urls(self, urls:List[str], max:int=0, format:str="", force:bool=False, workers:int=-1, tmp_dir:str|Path|None=None) -> None:
        """Downloads and adds sequences from a list of URLs to the SeqBank.

        This method processes a list of URLs, downloads the corresponding sequence files, and adds them to the SeqBank.
        It filters out URLs that have already been processed unless `force=True` is specified, and it can limit the number 
        of URLs processed based on the `max` argument. The processing can be parallelized using the `workers` argument.

        Args:
            urls (List[str]): A list of URLs to download and process.
            max (int, optional): Maximum number of URLs to process. If set to 0, all URLs will be processed. Defaults to 0.
            format (str, optional): The format of the sequence files (e.g., "fasta", "fastq"). If not provided, it will be auto-detected. Defaults to "".
            force (bool, optional): Whether to force re-processing of URLs even if they were processed before. Defaults to False.
            workers (int, optional): Number of workers to use for parallel processing. If set to -1, all available CPU cores will be used. Defaults to -1.
            tmp_dir (Union[str, Path, None], optional): A temporary directory to store downloaded files. Defaults to None.

        Returns:
            None
        """
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

    def ls(self) -> None:
        """
        Lists all accessions in the SeqBank.

        Iterates through the keys in the SeqBank and prints each one, decoded from bytes to ASCII.

        Returns:
            None
        """
        for k in self.file.keys():
            print(k.decode("ascii"))

    def add_files(self, files:List[str], max:int=0, format:str="", workers:int=1, filter:Path|list|set|None=None) -> None:
        """
        Adds sequences from multiple files to the SeqBank.

        This method processes a list of file paths, downloading and adding sequences from each file to the SeqBank.
        It supports parallel processing and optional filtering of specific accessions. The `max` argument limits the number of files to process.

        Args:
            files (List[str]): A list of file paths to process.
            max (int, optional): Maximum number of files to process. If set to 0, all files will be processed. Defaults to 0.
            format (str, optional): The format of the sequence files (e.g., "fasta", "fastq"). If not provided, it will be auto-detected. Defaults to "".
            workers (int, optional): Number of workers to use for parallel processing. Defaults to 1.
            filter (Union[Path, List[str], Set[str], None], optional): A filter for selecting specific accessions. Defaults to None.

        Returns:
            None
        """
        filter = parse_filter(filter)

        with Progress(*Progress.get_default_columns(), TimeElapsedColumn(), MofNCompleteColumn()) as progress:
            parallel = Parallel(n_jobs=workers, prefer="threads")
            add_file = delayed(self.add_file)
            overall_task = progress.add_task(f"[bold red]Adding files", total=len(files))
            parallel(add_file(file, progress=progress, format=format, overall_task=overall_task, filter=filter) for file in files)

    def copy(self, other: 'SeqBank') -> None:
        """
        Copies all entries from the current SeqBank to another SeqBank instance.

        This method iterates over all key-value pairs in the current SeqBank and adds them to the `other` SeqBank instance.
        The `other` SeqBank must be writable.

        Args:
            other (SeqBank): The target SeqBank instance where entries will be copied.

        Returns:
            None
        """
        for k, v in track(self.file.items()):
            other.file[k] = v

    def numpy(self, accession: str) -> np.ndarray:
        """
        Retrieves the sequence data for a given accession and returns it as an unsigned char NumPy array.

        Args:
            accession (str): The accession key for which the sequence data is retrieved.

        Returns:
            np.ndarray: The sequence data associated with the given accession, represented as an unsigned char NumPy array.
        """
        return np.frombuffer(self[accession], dtype="u1")

    def string(self, accession: str) -> str:
        """
        Retrieves the sequence data for a given accession and returns it as a string.

        Args:
            accession (str): The accession key for which the sequence data is retrieved.

        Returns:
            str: The sequence data associated with the given accession, represented as a string.
        """
        data = self[accession]
        return bytes_to_str(data)

    def record(self, accession: str) -> SeqRecord:
        """
        Retrieves the sequence data for a given accession and returns it as a BioPython SeqRecord object.

        Args:
            accession (str): The accession key for which the sequence data is retrieved.

        Returns:
            SeqRecord: A BioPython SeqRecord object containing the sequence data, with the given accession as its ID and an empty description.
        """
        record = SeqRecord(
            Seq(self.string(accession)),
            id=accession,
            description="",
        )
        return record

    def export(self, output: Path | str, format: str = "", accessions: List[str] | str | Path | None = None) -> None:
        """
        Exports the data from the SeqBank to a file using BioPython's SeqIO.

        Args:
            output (Path | str): The path or filename where the data should be exported.
            format (str, optional): The file format for exporting. If not specified, it will be inferred from the file extension.
            accessions (List[str] | str | Path | None, optional): A list of accessions to export. If a file path or string is provided, it will be read to obtain the list of accessions. If None, all accessions in the SeqBank are exported.

        Returns:
            None
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

        Returns:
            dict[str, int]: A dictionary mapping each accession to the length of its corresponding sequence.
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

        Args:
            nbins (int): The number of bins for the histogram. Default is 30.

        Returns:
            go.Figure: A Plotly figure object representing the histogram of sequence lengths.
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