"""Reproducible SeqBank engine benchmark; run with an environment containing all extras.

    python benchmarks/storage.py --output benchmarks/results.json

Uses real SeqBank FASTA imports, public reads, and direct adapter single-record
writes. Each engine/trial runs in a fresh process, sequentially. Datasets are
synthetic DNA; reads are warm-cache, not disk-cold. Durability defaults differ:
LMDB syncs each transaction; RocksDict/Speedict use default non-sync WAL writes.
"""
import argparse
import contextlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
import platform
import random
import resource
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

ENGINES = ('speedict', 'rocksdict', 'lmdb')
WORKLOADS = ((20000, 1000), (10000, 10000), (1000, 100000))


def worker(args):
    from seqbank import SeqBank
    bank_path = Path(args.workdir) / 'bank'
    bank = SeqBank(bank_path, write=True, engine=args.engine)
    # Include opening/initialization and closing/flushing in import time.
    start = time.perf_counter()
    with contextlib.redirect_stdout(io.StringIO()):
        bank.add_file(Path(args.fasta))
    bank.close()
    import_s = time.perf_counter() - start
    files = [p for p in bank_path.rglob('*') if p.is_file()]
    logical_size = sum(p.stat().st_size for p in files)
    allocated_size = sum(p.stat().st_blocks * 512 for p in files)
    reader = SeqBank(bank_path)
    keys = [f'seq{i:08d}' for i in range(args.count)]
    # Warm all records, verify contents, and exclude initial open from timing.
    checksum = 0
    for key in keys:
        value = reader[key]
        assert len(value) == args.length
        checksum += value[0]
    rng = random.Random(42)
    lookups = [keys[rng.randrange(len(keys))] for _ in range(20000)]
    start = time.perf_counter()
    for key in lookups:
        value = reader[key]
        checksum += value[0]
    lookup_s = time.perf_counter() - start
    start = time.perf_counter()
    total = sum(len(value) for _, value in reader.file.items())
    scan_s = time.perf_counter() - start
    assert total == args.count * args.length
    reader.close()
    # Expose the cost of the public add() API outside a batch transaction.
    single = SeqBank(Path(args.workdir) / 'single', write=True, engine=args.engine)
    single.file
    value = 'ACGT' * 250
    start = time.perf_counter()
    for i in range(1000):
        single.add(value, str(i))
    single.close()
    single_s = time.perf_counter() - start
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != 'darwin':
        rss *= 1024
    print(json.dumps(dict(
        engine=args.engine, records=args.count, sequence_length=args.length,
        import_seconds=import_s, import_mb_s=total / import_s / 1e6,
        random_reads_s=20000 / lookup_s, scan_mb_s=total / scan_s / 1e6,
        single_writes_s=1000 / single_s, disk_bytes=logical_size,
        allocated_bytes=allocated_size, peak_rss_bytes=rss, checksum=checksum,
    )))


def main(args):
    results = []
    with tempfile.TemporaryDirectory(prefix='seqbank-benchmark-') as temporary:
        root = Path(temporary)
        rng = random.Random(12345)
        for count, length in WORKLOADS:
            fasta = root / 'input.fasta'
            with fasta.open('w') as output:
                for i in range(count):
                    sequence = rng.choices('ACGTN', weights=(24, 24, 24, 24, 4), k=length)
                    output.write(f'>seq{i:08d}\n{"".join(sequence)}\n')
            for trial in range(args.repeats):
                # Rotate engine order to reduce ordering/temperature bias.
                order = ENGINES[trial % 3:] + ENGINES[:trial % 3]
                for engine in order:
                    directory = root / f'{engine}-{trial}'
                    directory.mkdir()
                    command = [sys.executable, str(Path(__file__).resolve()), '--worker',
                               '--engine', engine, '--workdir', str(directory),
                               '--fasta', str(fasta), '--count', str(count), '--length', str(length)]
                    completed = subprocess.run(command, capture_output=True, text=True, check=True)
                    result = json.loads(completed.stdout.strip().splitlines()[-1])
                    result['trial'] = trial + 1
                    results.append(result)
                    print(f'{count} x {length}, {engine}, trial {trial+1}: '
                          f'{result["import_mb_s"]:.1f} MB/s import, '
                          f'{result["random_reads_s"]:,.0f} reads/s', flush=True)
                    shutil.rmtree(directory)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata = dict(platform=platform.platform(), machine=platform.machine(),
                    python=sys.version, cpu_count=os.cpu_count(), repeats=args.repeats,
                    versions={engine: importlib.metadata.version(engine) for engine in ENGINES},
                    notes='Synthetic DNA; warm-cache reads; default engine durability; '
                          'FASTA import includes close; no forced compaction; MB = 1,000,000 bytes.')
    output.write_text(json.dumps(dict(environment=metadata, trials=results), indent=2) + '\n')
    lines = ['# SeqBank storage benchmark', '',
             f'Python {platform.python_version()}, {platform.machine()}, {args.repeats} trials per cell; medians.', '',
             '| Dataset | Engine | Import MB/s | Random reads/s | Scan MB/s | Single writes/s | Disk MB | Peak RSS MB |',
             '|---|---|---:|---:|---:|---:|---:|---:|']
    for count, length in WORKLOADS:
        for engine in ENGINES:
            subset = [r for r in results if r['records'] == count and r['engine'] == engine]
            median = lambda name: statistics.median(r[name] for r in subset)
            lines.append(f'| {count:,} × {length:,} bases | {engine} | '
                         f'{median("import_mb_s"):.1f} | {median("random_reads_s"):,.0f} | '
                         f'{median("scan_mb_s"):.1f} | {median("single_writes_s"):,.0f} | '
                         f'{median("disk_bytes")/1e6:.1f} | {median("peak_rss_bytes")/1e6:.1f} |')
    lines += ['', metadata['notes'], '',
              'Imports use SeqBank’s current batches (at most 1,000 records or about 8 MiB).',
              'Single writes are 1,000 individual 1 KB add() calls; these have unequal durability:',
              'LMDB commits synchronously, while RocksDict/Speedict use their default non-sync WAL.',
              'Do not interpret that column as an equal-durability engine comparison.',
              'Disk size is after close, before reads; no forced compaction. RSS includes Python,',
              'SeqBank imports, caches, and resident memory mappings; it is not engine-only heap usage.',
              'These datasets fit in memory; results do not establish cold-read or huge-bank performance.', '']
    output.with_suffix('.md').write_text('\n'.join(lines))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='benchmarks/results.json')
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--engine', choices=ENGINES)
    parser.add_argument('--workdir')
    parser.add_argument('--fasta')
    parser.add_argument('--count', type=int)
    parser.add_argument('--length', type=int)
    arguments = parser.parse_args()
    worker(arguments) if arguments.worker else main(arguments)
