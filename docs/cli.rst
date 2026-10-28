Command Line Interface Reference
================================

Use ``seqbank --help`` to list commands and ``seqbank COMMAND --help`` for
complete arguments and options for the installed version.

Engine selection
----------------

``add``, ``add-sequence-from-file``, ``url``, ``refseq``, and ``dfam`` accept
``--engine rocksdict``, ``--engine lmdb``, or ``--engine speedict`` when creating
a bank. If omitted, a new bank uses RocksDict if it is installed, otherwise LMDB. Existing banks are detected
automatically; an explicit engine must match their recorded engine.

LMDB is always installed. RocksDict requires the ``seqbank[rocksdict]`` extra and
Speedict the ``seqbank[speedict]`` extra (Python 3.10–3.12); see :ref:`storage-engines`.

.. code-block:: bash

    seqbank add sequences.sb sequences.fasta --engine lmdb
    seqbank add sequences.sb additional.fasta
    seqbank url sequences.sb https://example.com/sequences.fasta --workers 4
    seqbank count sequences.sb

Conversion and copying
----------------------

``convert SOURCE DESTINATION --engine ENGINE`` copies every record, including
internal URL tracking records, into the chosen engine. The destination must not exist; Speedict conversion requires its optional
package. The source is preserved.
If conversion fails, remove the incomplete destination before retrying.

.. code-block:: bash

    seqbank convert old.sb new.sb --engine lmdb

``merge SOURCE DESTINATION`` merges all records into an existing bank using its recorded engine,
or creates a new bank (RocksDict if installed, otherwise LMDB). Use ``convert`` to choose a destination engine.
Neither command performs a byte-for-byte filesystem backup.

Commands
--------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Command
     - Purpose
   * - ``add PATH FILES...``
     - Import sequence files; accepts ``--format`` and ``--filter``.
   * - ``add-sequence-from-file PATH ACCESSION FILE``
     - Import a file containing one sequence under the chosen accession.
   * - ``url PATH URLS...``
     - Download and import sequence files; accepts ``--workers`` and ``--tmp-dir``.
   * - ``refseq PATH``
     - Download RefSeq files; accepts ``--workers``, ``--max``, and ``--tmp-dir``.
   * - ``dfam PATH``
     - Download DFam; accepts ``--release`` and ``--curated/--no-curated``.
   * - ``ls PATH``
     - List accessions.
   * - ``count PATH``
     - Count stored records, including internal tracking records.
   * - ``delete PATH ACCESSIONS...``
     - Delete the specified accessions.
   * - ``merge SOURCE DESTINATION``
     - Merge records from one bank into another.
   * - ``convert SOURCE DESTINATION --engine ENGINE``
     - Convert to a new bank using a selected engine.
   * - ``export PATH OUTPUT``
     - Export sequences; accepts ``--format`` (default: FASTA).
   * - ``histogram PATH``
     - Plot sequence lengths; accepts ``--output-path``, ``--show``, ``--nbins``,
       ``--min``, and ``--max``. Static image output requires Chrome or Chromium.
