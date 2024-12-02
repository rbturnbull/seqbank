import tempfile
from pathlib import Path
from seqbank import SeqBank, SeqBankError
import numpy as np
from Bio.SeqRecord import SeqRecord
import pytest
import plotly.graph_objs as go
from unittest.mock import patch, MagicMock


TEST_DATA_PATH = Path(__file__).parent / "testdata"


def test_seqbank_add():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)
        seqbank = SeqBank(path=tmpdirname/"seqbank.sb", write=True)
        seqbank.add("ATCG", "test")
        assert "test" in seqbank
        assert np.all(seqbank.numpy("test") == np.array([1, 4, 2, 3], dtype="u1"))


def test_seqbank_add_file():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)
        seqbank = SeqBank(path=tmpdirname/"seqbank.sb", write=True)

        seqbank.add_file(TEST_DATA_PATH/"NC_024664.1.trunc.fasta")
        assert "NC_024664.1" in seqbank
        assert len(seqbank["NC_024664.1"]) == 840

        seqbank.add_file(TEST_DATA_PATH/"NZ_JAJNFP010000161.1.fasta")
        assert len(seqbank["NZ_JAJNFP010000161.1"]) == 849

        seqbank.add_file(TEST_DATA_PATH/"NC_036112.1.fasta")
        assert len(seqbank["NC_036112.1"]) == 980
        
        seqbank.add_file(TEST_DATA_PATH/"NC_044840.1.fasta")
        assert len(seqbank["NC_044840.1"]) == 700

        seqbank.add_file(TEST_DATA_PATH/"NC_010663.1.gb", format="genbank")
        assert len(seqbank["NC_010663.1"]) == 1066

        seqbank.add_file(TEST_DATA_PATH/"NC_036113.1.fasta")
        assert len(seqbank["NC_036113.1"]) == 1050


def test_add_sequence_from_file():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)
        seqbank = SeqBank(path=tmpdirname/"seqbank.sb", write=True)

        seqbank.add_sequence_from_file("NC_024664", TEST_DATA_PATH/"NC_024664.1.trunc.fasta")
        assert "NC_024664" in seqbank
        assert "NC_024664.1" not in seqbank
        assert len(seqbank["NC_024664"]) == 840

        seqbank.add_sequence_from_file("NC_010663", TEST_DATA_PATH/"NC_010663.1.gb")
        assert "NC_010663.1" not in seqbank
        assert len(seqbank["NC_010663"]) == 1066

        # check assert fails for multiple sequences
        with pytest.raises(SeqBankError):
            seqbank.add_sequence_from_file("5seq", TEST_DATA_PATH/"5seq.fa.gz")


def test_seqbank_add_file_fasta_with_filter_skip():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)
        seqbank = SeqBank(path=tmpdirname/"seqbank.sb", write=True)
        
        # Create a filter that excludes 'NC_024664.1' to ensure it's skipped
        filter_set = {"NC_036113.1"}  # 'NC_036113.1' is in the filter, but 'NC_024664.1' is not
        
        # Add a fasta file
        seqbank.add_file(TEST_DATA_PATH/"NC_024664.1.trunc.fasta", filter=filter_set)
        
        # Assert that 'NC_024664.1' was skipped (since it's not in the filter)
        assert "NC_024664.1" not in seqbank
        
        # Add another file with an accession in the filter
        seqbank.add_file(TEST_DATA_PATH/"NC_036113.1.fasta", filter=filter_set)
        
        # Assert that 'NC_036113.1' was added
        assert "NC_036113.1" in seqbank
        assert len(seqbank["NC_036113.1"]) == 1050

def test_seqbank_add_file_genbank_with_filter_skip():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)
        seqbank = SeqBank(path=tmpdirname/"seqbank.sb", write=True)
        
        # Create a filter that excludes 'NC_010663.1' to ensure it's skipped
        filter_set = {"NC_036112.1"}  # 'NC_036112.1' is in the filter, but 'NC_010663.1' is not
        
        # Add a genbank file
        seqbank.add_file(TEST_DATA_PATH/"NC_010663.1.gb", format="genbank", filter=filter_set)
        
        # Assert that 'NC_010663.1' was skipped (since it's not in the filter)
        assert "NC_010663.1" not in seqbank
        
        # Add another file with a record ID in the filter
        seqbank.add_file(TEST_DATA_PATH/"NC_036112.1.fasta", filter=filter_set)
        
        # Assert that 'NC_036112.1' was added
        assert "NC_036112.1" in seqbank
        assert len(seqbank["NC_036112.1"]) == 980

        
def test_record():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)
        seqbank = SeqBank(path=tmpdirname/"seqbank.sb", write=True)
        seqbank.add_file(TEST_DATA_PATH/"NC_036113.1.fasta")

        record = seqbank.record("NC_036113.1")
        

def test_string():
    seqbank = SeqBank(path=TEST_DATA_PATH/"seqbank.sb")
    
    string = seqbank.string("NC_036113.1")
    assert isinstance(string, str)
    assert len(string) == 1050
    assert string.startswith('GAAATACCCAATATCTTGTTCCAACAAGATAT')
    assert string.endswith('AATAAGGTAGGGATCATTAACACACC')
        

def test_record():
    seqbank = SeqBank(path=TEST_DATA_PATH/"seqbank.sb")

    record = seqbank.record("NC_010663.1")
    assert len(record) == 1066
    assert len(record.seq) == 1066
    assert isinstance(record, SeqRecord)

    assert record.id == "NC_010663.1"


def test_export():
    seqbank = SeqBank(path=TEST_DATA_PATH/"seqbank.sb")

    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)

        exported_path = tmpdirname/"exported.fasta"
        seqbank.export(exported_path)

        assert exported_path.exists()
        text = exported_path.read_text()
        assert text.count(">") == 6
        assert ">NC_024664.1\nTAAAAAGAAAAAAATTTT" in text
        assert ">NZ_JAJNFP010000161.1\nTAATATTTGCTTTCATTTCTAAATAG" in text
        assert ">NC_010663.1\nAGTTTTAAACCTCTGATCGAAC" in text
        assert ">NC_036112.1\nGAAATACCCAATATC" in text
        assert ">NC_036113.1\nGAAATACCCAATATCTTGTTCC" in text
        assert ">NC_044840.1\nGGGCGAACGACGGGAATTGAACCCGC" in text
        assert len(text) > 5000


def test_export_tsv_with_accessions():
    seqbank = SeqBank(path=TEST_DATA_PATH/"seqbank.sb")

    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)

        exported_path_tsv = tmpdirname/"exported.tsv"
        seqbank.export(exported_path_tsv, accessions=["NC_024664.1", "NC_010663.1"])

        assert exported_path_tsv.exists()
        text = exported_path_tsv.read_text()
        assert text.count("\t") == 2
        assert 'NC_024664.1\tTAAAAAGAAAAA' in text
        assert 'NC_010663.1\tAGTTTTAAAC' in text


@pytest.fixture
def prepare_accession_file():
    # Create a temporary file with a list of accessions
    with tempfile.NamedTemporaryFile(delete=False, mode='w') as f:
        f.write("NC_024664.1\nNZ_JAJNFP010000161.1\nNC_010663.1\n")
        accession_file = Path(f.name)
    yield accession_file
    # Clean up
    accession_file.unlink()


def test_export_from_file(prepare_accession_file):
    seqbank = SeqBank(path=TEST_DATA_PATH/"seqbank.sb")

    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)
        exported_path = tmpdirname/"exported.fasta"

        # Export using the file with accessions
        seqbank.export(exported_path, accessions=prepare_accession_file)

        # Validate the exported file
        assert exported_path.exists()
        text = exported_path.read_text()
        assert ">NC_024664.1\nTAAAAAGAAAAAAATTTT" in text
        assert ">NZ_JAJNFP010000161.1\nTAATATTTGCTTTCATTTCTAAATAG" in text
        assert ">NC_010663.1\nAGTTTTAAACCTCTGATCGAAC" in text
        assert len(text) > 100  # Adjust the length check based on actual content


def test_get_acccessions():
    seqbank = SeqBank(path=TEST_DATA_PATH/"seqbank.sb")
    accessions = seqbank.get_accessions()
    assert accessions == {'NC_036112.1', 'NC_024664.1', 'NC_010663.1', 'NC_036113.1', 'NC_044840.1', 'NZ_JAJNFP010000161.1'}


def test_close():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)
        seqbank_path = tmpdirname / "seqbank.sb"
        
        # Initialize SeqBank instance
        seqbank = SeqBank(path=seqbank_path, write=True)
        
        # Add some data to ensure the file is being used
        seqbank.add("ATCG", "test")
        
        # Close the SeqBank instance
        seqbank.close()
        
        # Verify that no exceptions are raised during the close
        try:
            seqbank.close()
        except Exception:
            assert False, "Exception raised during close when it should not be"
        

def test_file_not_found_error():
    # Test when the file does not exist
    temp_dir = Path('non_existing_directory')
    temp_dir.mkdir(parents=True, exist_ok=True)  # Create temp_dir

    # File path that does not exist
    non_existing_path = temp_dir / 'non_existing_file.sb'
    
    with pytest.raises(FileNotFoundError, match=f"Cannot find SeqBank file at path: {non_existing_path}"):
        # Initialize SeqBank with write=False and a non-existing path
        SeqBank(path=non_existing_path, write=False)

def test_no_exception_when_file_exists():
    # Test when the file exists
    temp_file = Path('existing_file.sb')
    
    try:
        # Create a temporary file
        temp_file.touch()
        
        # Initialize SeqBank with write=False and an existing path
        seqbank = SeqBank(path=temp_file, write=False)
        
        # No exception should be raised
        assert seqbank.path == temp_file.expanduser()
    finally:
        # Clean up
        if temp_file.exists():
            temp_file.unlink()

class MockSeqBank(SeqBank):
    def __init__(self, path: Path, write: bool = False):
        # Override to avoid file existence check
        self.path = path
        self.write = write
        # Mock file content based on the path
        if 'nonexistent' in str(path):
            self.file = None
        else:
            self.file = {'valid_key': np.array([0, 1, 2, 3], dtype="u1")}

    def __attrs_post_init__(self):
        # Skip file existence check
        pass

    def key(self, accession: str) -> str:
        return accession

def test_get_item_successful_retrieval():
   # Initialize MockSeqBank with a mock file
    mock_seqbank = MockSeqBank(path=Path('mock.sb'), write=False)
    
    # Retrieve data using a valid key
    result = mock_seqbank['valid_key']
    
    # Expected data as a NumPy array
    expected_result = np.array([0, 1, 2, 3], dtype="u1")
    
    # Assert that the result matches the expected data
    assert np.array_equal(result, expected_result)

def test_get_item_raises_exception():
    # Initialize MockSeqBank with a mock file (that doesn't contain the invalid key)
    mock_seqbank = MockSeqBank(path=Path('mock.sb'), write=False)
    
    # Attempt to retrieve data using an invalid key (which should raise an exception)
    invalid_key = 'invalid_key'
    
    with pytest.raises(SeqBankError, match=f"Failed to read {invalid_key} in SeqBank {mock_seqbank.path}"):
        result = mock_seqbank[invalid_key]  # This should trigger an exception

def test_contains_handles_exception():
    # Initialize MockSeqBank with a mock file
    mock_seqbank = MockSeqBank(path=Path('mock.sb'), write=False)
    
    # Mock an invalid scenario (e.g., accessing a nonexistent file)
    mock_seqbank.file = None  # Simulate a case where self.file is invalid
    
    # Use an arbitrary accession key, since it should trigger an exception due to mock_seqbank.file being None
    accession = "invalid_accession"
    
    # Call __contains__ and expect it to return False when an exception occurs
    result = accession in mock_seqbank
    
    # Assert that the result is False as expected
    assert result is False

# Mock SeqBank for testing
class MockSeqBankForItems(SeqBank):
    def __init__(self, path: Path, write: bool = False):
        super().__init__(path, write)
        # Mock file with sample data
        self.file = {
            'key1': b'\x01\x02\x03',
            'key2': b'\x04\x05\x06'
        }

    def __attrs_post_init__(self):
        # Skip file existence check
        pass

def test_items():
    mock_seqbank = MockSeqBankForItems(path=Path('mock.sb'), write=False)
    
    expected_items = {
        'key1': np.array([1, 2, 3], dtype='u1'),
        'key2': np.array([4, 5, 6], dtype='u1')
    }
    
    # Convert generator to dict
    result = dict(mock_seqbank.items())
    
    for key, expected_array in expected_items.items():
        assert key in result
        assert np.array_equal(result[key], expected_array)

def test_delete_key_error():
    # Initialize MockSeqBank with a mock file
    mock_seqbank = MockSeqBank(path=Path('mock.sb'), write=False)
    
    # Mock the file to simulate it contains only valid_key
    mock_seqbank.file = {'valid_key': np.array([0, 1, 2, 3], dtype="u1")}
    
    # Try deleting a key that does not exist in the mock file (this should raise a KeyError)
    invalid_key = 'invalid_key'
    
    # Call delete and ensure it doesn't raise an exception (even with KeyError inside the method)
    mock_seqbank.delete(invalid_key)

@pytest.fixture
def seqbank_with_data(tmp_path):
    """Fixture to set up a SeqBank instance with some predefined data."""
    seqbank_path = tmp_path / "seqbank.db"
    seqbank = SeqBank(path=seqbank_path, write=True)
    
    # Adding some sequences to SeqBank
    seqbank.add("ATCG", "seq1")
    seqbank.add("GCTA", "seq2")
    
    return seqbank   

def test_missing(seqbank_with_data):
    seqbank = seqbank_with_data
    
    # List of accessions to check
    accessions_to_check = ["seq1", "seq2", "seq3"]
    
    # Expected missing accessions
    expected_missing = {"seq3"}
    
    # Call the missing method without fetching missing accessions
    missing_accessions = seqbank.missing(accessions_to_check)
    
    assert missing_accessions == expected_missing


@pytest.fixture
def setup_seqbank(tmp_path):
    # Create a temporary directory for the SeqBank
    seqbank_path = tmp_path / "test_seqbank"
    seqbank_path.mkdir()

    # Create a SeqBank instance
    seqbank = SeqBank(path=seqbank_path, write=True)
    
    # Add sample sequences to the SeqBank
    sample_sequences = {
        "seq1": "ATCG",
        "seq2": "ATCGATCG",
        "seq3": "ATCGATCGATCG"
    }

    for accession, sequence in sample_sequences.items():
        seqbank.add(sequence, accession)

    return seqbank

def test_lengths_dict(setup_seqbank):
    # Retrieve the SeqBank instance from the fixture
    seqbank = setup_seqbank

    # Expected lengths of the sequences
    expected_lengths = {
        "seq1": 4,
        "seq2": 8,
        "seq3": 12
    }

    # Get lengths from the SeqBank
    lengths = seqbank.lengths_dict()

    # Assert that the lengths match the expected values
    assert lengths == expected_lengths

def test_histogram(setup_seqbank):
    # Retrieve the SeqBank instance from the fixture
    seqbank = setup_seqbank

    # Generate the histogram figure
    fig = seqbank.histogram()

    # Check that fig is a Plotly Figure object
    assert isinstance(fig, go.Figure)

    # Check that there is a trace in the figure with type 'histogram'
    histogram_traces = [trace for trace in fig.data if trace.type == 'histogram']

    # Assert that there is at least one histogram trace
    assert len(histogram_traces) == 1, "The figure does not contain one histogram trace."

# Test setup
@pytest.fixture
def seqbank_for_add():
    # Create a mock SeqBank instance
    with tempfile.TemporaryDirectory() as tmpdirname:
        seqbank = SeqBank(path=Path(tmpdirname) / "seqbank_mock.sb", write=True)
        # Prepopulate the SeqBank with some accessions
        seqbank.file = {'acc1': np.array([0, 1, 2, 3], dtype="u1")}
        return seqbank


def test_add_url_existing_url_no_force():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)
        seqbank = SeqBank(path=tmpdirname/"seqbank.sb", write=True)
        
        # Mocking the file to include the URL
        seqbank.file = {seqbank.key_url("http://example.com/file.fasta"): np.array([1, 2, 3], dtype="u1")}
        
        # Call add_url with the URL already present in the file and force=False
        result = seqbank.add_url("http://example.com/file.fasta", force=False)
        
        # Assert that the URL was not processed again
        assert not result

def test_add_url_exception_handling():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)
        seqbank = SeqBank(path=tmpdirname/"seqbank.sb", write=True)

        # Use patch to mock download_file to raise an exception
        with patch('seqbank.seqbank.download_file') as mock_download_file:
            mock_download_file.side_effect = Exception("Download failed")

            result = seqbank.add_url("http://example.com/file.fasta", tmp_dir=tmpdirname)
            
            # Assert that the exception was handled and False was returned
            assert not result
            
            # Verify the download_file was called with the expected arguments
            expected_local_path = tmpdirname / Path("file.fasta")
            # Get the actual call arguments
            actual_call_args = mock_download_file.call_args[0]

            # Extract and assert the URL and local path
            assert actual_call_args[0] == "http://example.com/file.fasta"
            assert actual_call_args[1].name == "file.fasta"  # Check if filename matches

def test_missing_exception_handling():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdirname = Path(tmpdirname)
        seqbank = SeqBank(path=tmpdirname / "seqbank.sb", write=True)

        # Define the mock behavior
        def mock_getitem(accession):
            if accession == "faulty_accession":
                raise ValueError("Data retrieval error")
            elif accession == "valid_accession":
                return np.array([0, 1, 2, 3], dtype="u1")
            # Simulate the case where accessions are not found
            raise KeyError("Accession not found")
        
        with patch.object(seqbank, '__getitem__', side_effect=mock_getitem):
            accessions = ["valid_accession", "faulty_accession", "missing_accession"]
            missing_accessions = seqbank.missing(accessions)
            
            # Check that 'faulty_accession' is in missing due to exception
            assert "faulty_accession" in missing_accessions
            
            # Check that 'missing_accession' is in missing (not present in mock)
            assert "missing_accession" in missing_accessions

@pytest.fixture
def seqbank():
    # Create a mock SeqBank instance
    with tempfile.TemporaryDirectory() as tmpdirname:
        seqbank = SeqBank(path=Path(tmpdirname) / "seqbank.sb", write=True)
        return seqbank

@patch.object(SeqBank, 'add_url')
@patch('joblib.Parallel')
@patch('joblib.delayed')
def test_add_urls_max(mock_delayed, mock_parallel, mock_add_url, seqbank):
    # Prepare mock data
    urls = [f'http://example.com/file{i}.fasta' for i in range(10)]  # 10 URLs
    max_urls = 5  # Limit to 5 URLs
    format = "fasta"
    force = False
    workers = 2
    tmp_dir = None

    # Mock the delayed function and parallel execution
    mock_delayed.return_value = lambda x: x
    mock_parallel.return_value = []

    # Call add_urls
    seqbank.add_urls(urls, max=max_urls, format=format, force=force, workers=workers, tmp_dir=tmp_dir)

    # Verify add_url is called only for the first `max_urls` URLs
    assert mock_add_url.call_count == max_urls

    # Check the URLs that were passed to add_url
    called_urls = [call[0][0] for call in mock_add_url.call_args_list]
    assert called_urls == urls[:max_urls]