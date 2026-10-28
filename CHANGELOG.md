# Release notes

## 0.3.0

- RocksDict is now an optional extra: `pip install 'seqbank[rocksdict]'`. A default
  `pip install seqbank` includes LMDB only. The `all` extra installs RocksDict and Speedict.
- New banks use RocksDict if it is installed and fall back to LMDB otherwise, so existing
  environments keep their behaviour. Existing banks always open with the engine recorded
  in `seqbank.json`; if it is not installed the error explains which extra to install.
  Use `seqbank convert` to move a bank to LMDB for use in environments without RocksDict.

## 0.2.0

- Add RocksDict and LMDB engines, installed by default with `pip install seqbank`.
  New banks default to RocksDict. Select another engine with `--engine` or the
  Python `Engine` enum; existing banks select their engine from `seqbank.json`.
- Retain legacy Speedict banks. Install `seqbank[speedict]` on Python 3.10–3.12
  to open or convert them. Speedict wheels are unavailable for Python 3.13–3.14.
- Add `seqbank convert SOURCE DESTINATION --engine lmdb` (or `rocksdict` or
  `speedict`). The destination must be new. Conversion preserves sequence and
  URL tracking records; it does not modify the source. Keep `seqbank.json` with
  the bank directory when copying or distributing banks.
- Batch file imports and conversion at 1,000 records or approximately 8 MiB.
  Individual `add()` calls and DFam HDF5 imports still commit per sequence.
- Support concurrent local LMDB writer processes with serialized transactions.
  Built-in parallel imports still use threads. RocksDict/Speedict require one
  writer process. Open database handles inside their owning process.
- Hash long URL tracking keys to fit LMDB key limits, including during conversion
  from legacy banks. Existing short URL keys remain unchanged.
- Support Python 3.10–3.14, with CLI engine validation and parameter help.
- Update Plotly/Kaleido compatibility. Static image export now requires Chrome
  or Chromium; install a browser with `uv run plotly_get_chrome -y` if needed.
- Include reproducible engine benchmarks. Results are warm-cache measurements
  with current engine defaults, which provide different durability guarantees.
