"""Byte-key/byte-value storage adapters and persistent engine identification."""
import importlib
import json
import hashlib
from enum import Enum
from threading import RLock

from .exceptions import SeqBankError

class Engine(str, Enum):
    """Supported database engines (also accepted as strings in the Python API)."""

    SPEEDICT = "speedict"
    ROCKSDICT = "rocksdict"
    LMDB = "lmdb"

    def __str__(self):
        return self.value


ENGINES = tuple(engine.value for engine in Engine)
METADATA = "seqbank.json"


def _storage_key(key):
    """Bound URL tracking keys for LMDB while preserving the internal namespace."""
    if key.startswith(b"/seqbank/url/") and len(key) > 511:
        return b"/seqbank/url-sha256/" + hashlib.sha256(key).hexdigest().encode("ascii")
    return key


class RocksStorage:
    def __init__(self, module, path, write):
        options = module.Options(raw_mode=True)
        options.set_compression_type(module.DBCompressionType.none())
        options.set_optimize_filters_for_hits(True)
        options.optimize_for_point_lookup(1024)
        options.set_max_open_files(500)
        self.module = module
        self.db = module.Rdict(
            str(path), options=options,
            access_type=module.AccessType.read_write() if write else module.AccessType.read_only(),
        )

    def __getitem__(self, key):
        try:
            return self.db[key]
        except KeyError:
            normalized = _storage_key(key)
            if normalized == key:
                raise
            return self.db[normalized]

    def __setitem__(self, key, value):
        self.db[self._write_key(key)] = bytes(value)

    def __delitem__(self, key):
        del self.db[self._write_key(key)]

    def __contains__(self, key):
        if key in self.db:
            return True
        normalized = _storage_key(key)
        return normalized != key and normalized in self.db

    def items(self):
        return self.db.items()

    def keys(self):
        return self.db.keys()

    def _write_key(self, key):
        normalized = _storage_key(key)
        # Preserve legacy raw URL records when updating an existing bank.
        return key if normalized != key and key in self.db else normalized

    def put_many(self, records):
        batch = self.module.WriteBatch(raw_mode=True)
        for key, value in records:
            batch.put(self._write_key(key), bytes(value))
        self.db.write(batch)

    def close(self):
        self.db.close()


class LMDBStorage:
    def __init__(self, module, path, write):
        # On 64-bit systems this reserves virtual address space, not RAM.
        self.db = module.open(str(path), readonly=not write, create=write, map_size=1 << 40)
        self._write_lock = RLock()

    def __getitem__(self, key):
        with self.db.begin() as txn:
            value = txn.get(_storage_key(key))
            if value is None:
                raise KeyError(key)
            return value

    def __setitem__(self, key, value):
        with self._write_lock, self.db.begin(write=True) as txn:
            txn.put(_storage_key(key), bytes(value))

    def __delitem__(self, key):
        with self._write_lock, self.db.begin(write=True) as txn:
            if not txn.delete(_storage_key(key)):
                raise KeyError(key)

    def __contains__(self, key):
        with self.db.begin() as txn:
            return txn.get(_storage_key(key)) is not None

    def items(self):
        with self.db.begin() as txn:
            with txn.cursor() as cursor:
                yield from cursor

    def keys(self):
        for key, _ in self.items():
            yield key

    def put_many(self, records):
        with self._write_lock, self.db.begin(write=True) as txn:
            for key, value in records:
                txn.put(_storage_key(key), bytes(value))

    def close(self):
        self.db.close()


def open_storage(path, write, requested=None):
    """Identify existing banks; never try a different engine on an existing bank."""
    if requested is not None and requested not in ENGINES:
        raise SeqBankError(f"Unknown engine {requested!r}; choose from {', '.join(ENGINES)}")
    if requested is not None:
        requested = Engine(requested)
    metadata_path = path / METADATA
    existing = path.exists() and any(path.iterdir())
    if metadata_path.exists():
        try:
            metadata = json.loads(metadata_path.read_text())
            engine = Engine(metadata["engine"])
            if metadata["format_version"] != 1 or engine not in ENGINES:
                raise ValueError("unsupported format version or engine")
        except (ValueError, KeyError, TypeError) as err:
            raise SeqBankError(f"Invalid SeqBank metadata at {metadata_path}: {err}") from err
    elif existing:
        if not (path / "CURRENT").exists():
            raise SeqBankError(f"Cannot identify SeqBank engine at {path}")
        # All banks created before engine metadata were Speedict banks.
        engine = Engine.SPEEDICT
    else:
        engine = requested or Engine.ROCKSDICT
    if requested is not None and requested != engine:
        raise SeqBankError(f"Bank uses {engine}, not {requested}; convert it to a new directory instead")
    try:
        module = importlib.import_module(engine.value)
    except ImportError as err:
        hint = " Speedict wheels require Python <3.13." if engine == "speedict" else ""
        install = "seqbank[speedict]" if engine == Engine.SPEEDICT else "seqbank"
        raise SeqBankError(
            f"Engine {engine.value!r} is unavailable. Install with: pip install '{install}'.{hint}"
        ) from err
    if write:
        path.mkdir(parents=True, exist_ok=True)
    storage = LMDBStorage(module, path, write) if engine == "lmdb" else RocksStorage(module, path, write)
    if write and not metadata_path.exists():
        try:
            # Exclusive creation also permits independent LMDB writer processes.
            with metadata_path.open("x") as output:
                json.dump({"format_version": 1, "engine": engine}, output)
        except FileExistsError:
            pass
        except Exception:
            storage.close()
            raise
    return storage, engine
