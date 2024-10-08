from seqbank.utils import parse_filter
from pathlib import Path

TEST_DATA_PATH = Path(__file__).parent / "testdata"

def test_parse_filter():
    assert parse_filter(None) == None
    assert parse_filter(["a", "b", "c"]) == {"a", "b", "c"}
    assert parse_filter({"a", "b", "c"}) == {"a", "b", "c"}
    assert parse_filter(TEST_DATA_PATH/"accessions.txt") == {
        "NC_036112.1",
        "NC_024664.1",
        "NC_010663",
        "NC_036113.1",
        "NC_044840.1",
        "NZ_JAJNFP010000161.1",
    }