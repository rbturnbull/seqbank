"""Byte-key/byte-value storage adapters and persistent engine identification."""

import importlib
import json
import hashlib
from enum import Enum
from pathlib import Path
from threading import RLock
from types import ModuleType
from typing import Iterable, Iterator, Protocol

from .exceptions import SeqBankError


class Engine(str, Enum):
    """Supported database engines (also accepted as strings in the Python API)."""

    SPEEDICT = "speedict"
    ROCKSDICT = "rocksdict"
    LMDB = "lmdb"

    def __str__(self) -> str:
        return self.value


ENGINES = tuple(engine.value for engine in Engine)
METADATA = "seqbank.json"


def _storage_key(key: bytes) -> bytes:
    """Bound URL tracking keys for LMDB while preserving the internal namespace."""
    if key.startswith(b"/seqbank/url/") and len(key) > 511:
        return b"/seqbank/url-sha256/" + hashlib.sha256(key).hexdigest().encode("ascii")
    return key


class Storage(Protocol):
    """Byte-key/byte-value interface implemented by every engine adapter."""

    def __getitem__(self, key: bytes) -> bytes: ...

    def __setitem__(self, key: bytes, value: bytes) -> None: ...

    def __delitem__(self, key: bytes) -> None: ...

    def __contains__(self, key: bytes) -> bool: ...

    def items(self) -> Iterator[tuple[bytes, bytes]]: ...

    def keys(self) -> Iterator[bytes]: ...

    def put_many(self, records: Iterable[tuple[bytes, bytes]]) -> None: ...

    def close(self) -> None: ...


class RocksStorage:
    def __init__(self, module: ModuleType, path: Path, write: bool) -> None:
        options = module.Options(raw_mode=True)
        options.set_compression_type(module.DBCompressionType.none())
        options.set_optimize_filters_for_hits(True)
        options.optimize_for_point_lookup(1024)
        options.set_max_open_files(500)
        self.module = module
        self.db = module.Rdict(
            str(path),
            options=options,
            access_type=module.AccessType.read_write() if write else module.AccessType.read_only(),
        )

    def __getitem__(self, key: bytes) -> bytes:
        try:
            return self.db[key]
        except KeyError:
            normalized = _storage_key(key)
            if normalized == key:
                raise
            return self.db[normalized]

    def __setitem__(self, key: bytes, value: bytes) -> None:
        self.db[self._write_key(key)] = bytes(value)

    def __delitem__(self, key: bytes) -> None:
        del self.db[self._write_key(key)]

    def __contains__(self, key: bytes) -> bool:
        if key in self.db:
            return True
        normalized = _storage_key(key)
        return normalized != key and normalized in self.db

    def items(self) -> Iterator[tuple[bytes, bytes]]:
        return self.db.items()

    def keys(self) -> Iterator[bytes]:
        return self.db.keys()

    def _write_key(self, key: bytes) -> bytes:
        normalized = _storage_key(key)
        # Preserve legacy raw URL records when updating an existing bank.
        return key if normalized != key and key in self.db else normalized

    def put_many(self, records: Iterable[tuple[bytes, bytes]]) -> None:
        batch = self.module.WriteBatch(raw_mode=True)
        for key, value in records:
            batch.put(self._write_key(key), bytes(value))
        self.db.write(batch)

    def close(self) -> None:
        self.db.close()


class LMDBStorage:
    def __init__(self, module: ModuleType, path: Path, write: bool) -> None:
        # On 64-bit systems this reserves virtual address space, not RAM.
        # Read-only banks skip LMDB's lock file: readers then take no file locks, which can hang
        # on network filesystems when many processes open the same bank.
        # while it is being read this way.
        self.db = module.open(str(path), readonly=not write, create=write, map_size=1 << 40, lock=write)
        self._write_lock = RLock()

    def __getitem__(self, key: bytes) -> bytes:
        with self.db.begin() as txn:
            value = txn.get(_storage_key(key))
            if value is None:
                raise KeyError(key)
            return value

    def __setitem__(self, key: bytes, value: bytes) -> None:
        with self._write_lock, self.db.begin(write=True) as txn:
            txn.put(_storage_key(key), bytes(value))

    def __delitem__(self, key: bytes) -> None:
        with self._write_lock, self.db.begin(write=True) as txn:
            if not txn.delete(_storage_key(key)):
                raise KeyError(key)

    def __contains__(self, key: bytes) -> bool:
        with self.db.begin() as txn:
            return txn.get(_storage_key(key)) is not None

    def items(self) -> Iterator[tuple[bytes, bytes]]:
        with self.db.begin() as txn:
            with txn.cursor() as cursor:
                yield from cursor

    def keys(self) -> Iterator[bytes]:
        for key, _ in self.items():
            yield key

    def put_many(self, records: Iterable[tuple[bytes, bytes]]) -> None:
        with self._write_lock, self.db.begin(write=True) as txn:
            for key, value in records:
                txn.put(_storage_key(key), bytes(value))

    def close(self) -> None:
        self.db.close()


INSTALL_HINTS = {
    Engine.ROCKSDICT: "pip install 'seqbank[rocksdict]'",
    Engine.SPEEDICT: "pip install 'seqbank[speedict]' (Speedict wheels require Python <3.13)",
    Engine.LMDB: "pip install seqbank",
}


def default_engine() -> Engine:
    """Engine for new banks: RocksDict if it can be imported, otherwise LMDB."""
    try:
        importlib.import_module(Engine.ROCKSDICT.value)
    except ImportError:
        return Engine.LMDB
    return Engine.ROCKSDICT


def open_storage(path: Path, write: bool, requested: Engine | str | None = None) -> tuple[Storage, Engine]:
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
        engine = requested or default_engine()
    if requested is not None and requested != engine:
        raise SeqBankError(f"Bank uses {engine}, not {requested}; convert it to a new directory instead")
    try:
        module = importlib.import_module(engine.value)
    except ImportError as err:
        raise SeqBankError(f"Engine {engine.value!r} is unavailable. Install with: {INSTALL_HINTS[engine]}.") from err
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
