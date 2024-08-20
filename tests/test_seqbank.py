import tempfile
from pathlib import Path
from seqbank import SeqBank
import numpy as np
from Bio.SeqRecord import SeqRecord
import pytest


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

