from seqbank.io import get_file_format, seq_count
import pytest
from pathlib import Path

TEST_DATA_PATH = Path(__file__).parent / "testdata"


def test_get_file_format():
    assert get_file_format("test.fasta") == "fasta"
    assert get_file_format("test.fa") == "fasta"
    assert get_file_format("test.fna") == "fasta"    
    assert get_file_format("test.fna.gz") == "fasta"        

    assert get_file_format("test.genbank") == "genbank"        
    assert get_file_format("test.genbank.gz") == "genbank"        
    assert get_file_format("test.gb.gz") == "genbank"        
    assert get_file_format("test.gbk.gz") == "genbank"        

    assert get_file_format("test.embl") == "embl"
    assert get_file_format("test.embl.gz") == "embl"
    assert get_file_format("test.embl.bz2") == "embl"

    assert get_file_format("test.nexus") == "nexus"
    assert get_file_format("test.nxs.gz") == "nexus"

    assert get_file_format("test.tsv.gz") == "tab"        
    
    assert get_file_format("test.fastq") == "fastq"        
    assert get_file_format("test.fq") == "fastq"        

    with pytest.raises(ValueError):
        get_file_format("test.txt")


def test_seq_count():
    seq_count(TEST_DATA_PATH/"NC_010663.1.gb") == 1
    seq_count(TEST_DATA_PATH/"NC_010663.1.gb.bz2") == 1

    seq_count(TEST_DATA_PATH/"5seq.fa.gz") == 5
    seq_count(TEST_DATA_PATH/"5seq.fa") == 5
    