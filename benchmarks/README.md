Run the benchmark in Python 3.10–3.12 with all three engines installed:

```sh
uv sync --all-extras
uv run python benchmarks/storage.py --output benchmarks/results.json
```

The script generates deterministic synthetic DNA in three sizes, runs three
trials per engine in separate subprocesses, and deletes temporary databases.
It saves every measurement and package version to JSON and median results to
Markdown. Engines run sequentially with their order rotated between trials.

Measurements cover SeqBank FASTA import (including database close), warm random
lookups through `SeqBank.__getitem__`, byte-level scans, unbatched public `add()`
calls, database size, and whole-process peak resident memory. Imported values
and total payload sizes are checked. No existing user databases are changed.

These results compare the current adapters, including their settings. They are
not an equal-durability or universal engine ranking: LMDB commits synchronously,
while the Speedict and RocksDict defaults do not sync each write to disk.
Reads are warm-cache and datasets fit in memory. Disk size is measured after
close without forced compaction; repeat timing on representative production
data before choosing a backend. The included run was on a disk with roughly
4.8 GiB free, which may affect write measurements.
