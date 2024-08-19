# -*- coding: utf-8 -*-

import numpy as np

vocab_to_int = {"A": 1, "C": 2, "G": 3, "T": 4, "N": 0}
int_to_vocab = dict(zip(vocab_to_int.values(), vocab_to_int.keys()))

DELETE = bytes(range(256)).translate(None, b'ACGTNacgtn')
TABLE = bytearray(b'\0' * 256)
REVERSE_TABLE = bytearray(b'\0' * 256)

for char, value in vocab_to_int.items():
    TABLE[ord(char)] = value
    TABLE[ord(char.lower())] = value

    REVERSE_TABLE[value] = ord(char)


def seq_to_bytes(seq:str|bytes, delete:bool=True) -> bytes:
    """
    Convert a DNA sequence to bytes.

    Args:
        seq (str | bytes): _description_
        delete (bool, optional): _description_. Defaults to True.

    Returns:
        bytes: _description_
    """
    if isinstance(seq, str):
        seq = seq.encode("ascii")

    assert isinstance(seq, bytes), f"Expected bytes, got {type(seq)}"

    if delete:
        return seq.translate(TABLE, DELETE)
    return seq.translate(TABLE)


def seq_to_numpy(seq:str|bytes, delete:bool=True) -> np.ndarray:
    b = seq_to_bytes(seq, delete=delete)
    return np.frombuffer(b, dtype="u1")


def bytes_to_str(seq:bytes) -> str:
    return seq.translate(REVERSE_TABLE).decode("ascii")

