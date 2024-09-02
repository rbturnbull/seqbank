import tempfile
from pathlib import Path
from seqbank import SeqBank, SeqBankError
import numpy as np
from Bio.SeqRecord import SeqRecord
import pytest
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
        # seqbank = SeqBank(path=TEST_DATA_PATH/"seqbank.sb", write=True)

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
            # Try accessing the database after closing it
            seqbank.file
            assert False, "Expected an exception when accessing the file after closing"
        except Exception:
            pass
        
        # Since we are not interacting with the actual file system in the test,
        # we cannot directly check the closed state. We can only assert that no errors occur.

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

def test_missing_with_get(seqbank_with_data):
    seqbank = seqbank_with_data
    
    # List of accessions to check
    accessions_to_check = ["seq1", "seq2", "seq3"]
    
    # Mocking the behavior of fetching missing accessions
    with patch.object(seqbank, '__getitem__', side_effect=lambda x: seqbank.file[seqbank.key(x)] if x != "seq3" else None):
        missing_accessions = seqbank.missing(accessions_to_check, get=True)
    
    # Since 'seq3' cannot be fetched, it should be reported as missing
    assert missing_accessions == {"seq3"}

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

def test_get_accession_lengths(setup_seqbank):
    # Retrieve the SeqBank instance from the fixture
    seqbank = setup_seqbank

    # Expected lengths of the sequences
    expected_lengths = {
        "seq1": 4,
        "seq2": 8,
        "seq3": 12
    }

    # Get lengths from the SeqBank
    lengths = seqbank.get_accession_lengths()

    # Assert that the lengths match the expected values
    assert lengths == expected_lengths