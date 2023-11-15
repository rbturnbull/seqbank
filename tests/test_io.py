from seqbank.io import get_file_format
import pytest


def test_get_file_format():
    assert get_file_format("test.fasta") == "fasta"
    assert get_file_format("test.fa") == "fasta"
    assert get_file_format("test.fna") == "fasta"    
    assert get_file_format("test.fna.gz") == "fasta"        


    assert get_file_format("test.genbank") == "genbank"        
    assert get_file_format("test.genbank.gz") == "genbank"        
    assert get_file_format("test.gb.gz") == "genbank"        
    assert get_file_format("test.gbk.gz") == "genbank"        


    assert get_file_format("test.tsv.gz") == "tab"        
    
    assert get_file_format("test.fastq") == "fastq"        
    assert get_file_format("test.fq") == "fastq"        

    with pytest.raises(ValueError):
        get_file_format("test.txt")