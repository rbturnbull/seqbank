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


def seq_to_bytes(seq: str | bytes, delete: bool = True) -> bytes:
    """
    Convert a DNA sequence to a byte representation using a predefined mapping.

    Args:
        seq (str | bytes): The DNA sequence to convert. Can be a string or bytes.
        delete (bool, optional): If True, delete characters that are not in 'ACGTN' or 'acgtn'. Defaults to True.

    Returns:
        bytes: The byte-encoded version of the input sequence.
    """
    if isinstance(seq, str):
        seq = seq.encode("ascii")

    assert isinstance(seq, bytes), f"Expected bytes, got {type(seq)}"

    if delete:
        return seq.translate(TABLE, DELETE)
    return seq.translate(TABLE)


def seq_to_numpy(seq: str | bytes, delete: bool = True) -> np.ndarray:
    """
    Convert a DNA sequence to a NumPy array of unsigned 8-bit integers.

    Args:
        seq (str | bytes): The DNA sequence to convert, in string or byte format.
        delete (bool, optional): If True, remove any characters not in 'ACGTN'. Defaults to True.

    Returns:
        np.ndarray: A NumPy array representing the byte-encoded DNA sequence.
    """
    b = seq_to_bytes(seq, delete=delete)
    return np.frombuffer(b, dtype="u1")


def bytes_to_str(seq: bytes) -> str:
    """
    Convert byte-encoded DNA sequence back to a string.

    Args:
        seq (bytes): The byte-encoded DNA sequence.

    Returns:
        str: The original DNA sequence in string format.
    """
    return seq.translate(REVERSE_TABLE).decode("ascii")