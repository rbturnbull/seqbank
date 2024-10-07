================
Seqbank
================

.. start-badges

|pypi badge| |testing badge| |coverage badge| |docs badge| |black badge| |git3moji badge|

.. |pypi badge| image:: https://img.shields.io/pypi/v/seqbank
    :target: https://pypi.org/project/seqbank/

.. |testing badge| image:: https://github.com/rbturnbull/seqbank/actions/workflows/testing.yml/badge.svg
    :target: https://github.com/rbturnbull/seqbank/actions

.. |docs badge| image:: https://github.com/rbturnbull/seqbank/actions/workflows/docs.yml/badge.svg
    :target: https://rbturnbull.github.io/seqbank
    
.. |black badge| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    
.. |coverage badge| image:: https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/rbturnbull/b1625e7f45428007f0982543d9d346d0/raw/coverage-badge.json
    :target: https://rbturnbull.github.io/seqbank/coverage/

.. |git3moji badge| image:: https://img.shields.io/badge/git3moji-%E2%9A%A1%EF%B8%8F%F0%9F%90%9B%F0%9F%93%BA%F0%9F%91%AE%F0%9F%94%A4-fffad8.svg
    :target: https://robinpokorny.github.io/git3moji/
        
.. end-badges

.. start-quickstart

SeqBank is a powerful and flexible command-line application designed to simplify the management and processing of large DNA sequence datasets. Whether you're working with local sequence files, 
retrieving data from remote URLs, or integrating sequences from databases like RefSeq and DFam, SeqBank provides an efficient, user-friendly solution.

SeqBank allows users to quickly add, organize, and manipulate sequences using a structured, numerical format optimized for fast retrieval and analysis. 
It's especially useful for bioinformatics professionals who regularly handle vast amounts of genomic data.


Installation
============

To install the latest version from the repository, you can use this command:

.. code-block:: bash

    pip install seqbank

Or install directly from the GitHub repository:

.. code-block:: bash

    pip install git+https://github.com/rbturnbull/seqbank.git

.. note ::

    Soon seqbank will be able to be installed using conda.

Usage
===========
    
SeqBank provides a command-line interface (CLI) for managing DNA sequence data efficiently. Below are the main tools, along with examples of how to use them in practical workflows.

Adding Sequences
----------------

SeqBank allows you to import sequence data from files or URLs into the database. The system supports multiple sequence formats, providing flexibility in handling various datasets.

**Example:**

To add sequences from one or more local files:

.. code-block:: bash

    seqbank add /path/to/seqbank /path/to/sequence1.fasta /path/to/sequence2.fasta --format fasta

To add sequences from a list of URLs:

.. code-block:: bash

    seqbank url /path/to/seqbank https://example.com/sequence1.fasta https://example.com/sequence2.fasta --format fasta --workers 4

**Use case:**  
Suppose you have a new set of genome sequences in FASTA format stored locally or accessible via URLs. You can quickly import these sequences into your SeqBank database for centralized storage and further analysis.


Managing Databases
------------------

SeqBank provides commands to manage and query the sequences in your database. You can list, count, and delete sequences, allowing efficient database management.

**Example:**

To list all sequences in the database:

.. code-block:: bash

    seqbank ls /path/to/seqbank

To count the number of sequences stored:

.. code-block:: bash

    seqbank count /path/to/seqbank

To delete a specific sequence by accession number:

.. code-block:: bash

    seqbank delete /path/to/seqbank ABC123DEF456

**Use case:**  
If you're managing a growing sequence database, the `ls` command can help you track the sequences, while `delete` can be used to remove outdated or incorrect entries.


Exporting Sequences
-------------------

You can export your stored sequences to common formats like FASTA for easy sharing and use with other bioinformatics tools. This ensures compatibility with external platforms.

**Example:**

To export sequences in FASTA format to a specific output directory:

.. code-block:: bash

    seqbank export /path/to/seqbank /output/directory --format fasta

**Use case:**  
After storing a collection of curated sequences, you may need to export them in FASTA format for downstream analysis using tools like BLAST or multiple sequence alignment software.


Integration with RefSeq and DFam
--------------------------------

SeqBank integrates with popular genomic databases like RefSeq and DFam, allowing users to download and incorporate sequences from these sources.

**Example:**

To download and add RefSeq sequences with a maximum of 1000 sequences using 4 workers:

.. code-block:: bash

    seqbank refseq /path/to/seqbank --max 1000 --workers 4

To download and add DFam sequences from the current release with curated data:

.. code-block:: bash

    seqbank dfam /path/to/seqbank --release current --curated

**Use case:**  
If you are studying repetitive elements in a genome, you can easily integrate sequences from DFam into your SeqBank database for comprehensive analysis.


Visualization of Sequence Data
------------------------------

SeqBank includes built-in functionality for generating histograms of sequence lengths, providing a visual summary of the data.

**Example:**

To generate and save a histogram of sequence lengths:

.. code-block:: bash

    seqbank histogram /path/to/seqbank --output histogram.png --nbins 50

To generate and display the histogram interactively:

.. code-block:: bash

    seqbank histogram /path/to/seqbank --show --nbins 50

**Use case:**  
When working with a dataset of varying sequence lengths, generating a histogram can help visualize the distribution and detect outliers or inconsistencies in the data.


Copying Databases
-----------------

SeqBank allows you to copy sequences from one SeqBank database to another, facilitating data migration or backup processes.

**Example:**

To copy sequences from a source SeqBank to a destination SeqBank:

.. code-block:: bash

    seqbank cp /path/to/source_seqbank /path/to/destination_seqbank

**Use case:**  
For maintaining backups of your sequence database or migrating data to a new location, the `cp` command provides a straightforward method to duplicate your SeqBank data.


Filtering Sequences and Custom Workflows
----------------------------------------

SeqBank supports filtering sequences based on criteria such as sequence length or file format before adding them to the database. Additionally, multi-threaded downloading allows you to download and process sequences more efficiently.

**Example:**

To filter sequences longer than 1000 bp before adding them:

.. code-block:: bash

    seqbank add /path/to/seqbank /path/to/sequences.fasta --format fasta --filter /path/to/filter_file

To enable multi-threaded downloading when adding sequences from URLs:

.. code-block:: bash

    seqbank url /path/to/seqbank https://example.com/sequence1.fasta https://example.com/sequence2.fasta --format fasta --workers 4 --tmp-dir /path/to/tmp

**Use case:**  
In projects where only sequences longer than a specific threshold are required, the filtering feature ensures that only relevant sequences are stored. Multi-threaded downloading can be utilized when processing large datasets to save time.

.. end-quickstart


Credits
==================================

.. start-credits

* Robert Turnbull <robert.turnbull@unimelb.edu.au>
* Rafsan Al Mamun <rafsan7238@gmail.com>

.. end-credits

