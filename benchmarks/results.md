# SeqBank storage benchmark

Python 3.11.11, arm64, 3 trials per cell; medians.

| Dataset | Engine | Import MB/s | Random reads/s | Scan MB/s | Single writes/s | Disk MB | Peak RSS MB |
|---|---|---:|---:|---:|---:|---:|---:|
| 20,000 × 1,000 bases | speedict | 196.0 | 630,812 | 2324.1 | 78,891 | 20.6 | 119.8 |
| 20,000 × 1,000 bases | rocksdict | 196.8 | 637,669 | 2513.2 | 55,866 | 20.6 | 120.6 |
| 20,000 × 1,000 bases | lmdb | 82.7 | 815,274 | 6866.0 | 118 | 22.1 | 97.4 |
| 10,000 × 10,000 bases | speedict | 537.6 | 337,217 | 5127.0 | 102,125 | 100.7 | 226.2 |
| 10,000 × 10,000 bases | rocksdict | 438.3 | 382,859 | 5442.3 | 104,435 | 100.7 | 221.7 |
| 10,000 × 10,000 bases | lmdb | 279.6 | 704,455 | 28565.3 | 118 | 164.3 | 259.4 |
| 1,000 × 100,000 bases | speedict | 617.8 | 98,301 | 9064.2 | 96,445 | 100.1 | 221.7 |
| 1,000 × 100,000 bases | rocksdict | 426.6 | 104,148 | 9802.5 | 99,442 | 100.1 | 224.7 |
| 1,000 × 100,000 bases | lmdb | 330.7 | 325,899 | 47356.9 | 118 | 114.9 | 228.9 |

Synthetic DNA; warm-cache reads; default engine durability; FASTA import includes close; no forced compaction; MB = 1,000,000 bytes.

Imports use SeqBank’s current batches (at most 1,000 records or about 8 MiB).
Single writes are 1,000 individual 1 KB add() calls; these have unequal durability:
LMDB commits synchronously, while RocksDict/Speedict use their default non-sync WAL.
Do not interpret that column as an equal-durability engine comparison.
Disk size is after close, before reads; no forced compaction. RSS includes Python,
SeqBank imports, caches, and resident memory mappings; it is not engine-only heap usage.
These datasets fit in memory; results do not establish cold-read or huge-bank performance.
