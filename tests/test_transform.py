from seqbank import transform
import numpy as np


def test_seq_to_bytes():
    assert transform.seq_to_bytes("ACGTN") == b'\x01\x02\x03\x04\x00'


def test_seq_to_numpy():
    assert (transform.seq_to_numpy("ACGTN") == np.array([1, 2, 3, 4, 0], dtype="u1")).all()


def test_seq_to_numpy_extra():
    assert (transform.seq_to_numpy("ACGTNJIKLO") == np.array([1, 2, 3, 4, 0], dtype="u1")).all()


def test_seq_to_numpy_extra_no_delete():
    assert (transform.seq_to_numpy("ACGTNJIKLO", delete=False) == np.array([1, 2, 3, 4, 0, 0, 0, 0, 0, 0], dtype="u1")).all()


def test_bytes_to_str():
    assert transform.bytes_to_str(b'\x01\x02\x03\x04\x00') == "ACGTN"
    assert transform.bytes_to_str(transform.seq_to_bytes("ACJIKLOGT")) == "ACGT"
