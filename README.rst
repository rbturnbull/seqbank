================
Seqbank
================

.. start-badges

|testing badge| |coverage badge| |docs badge| |black badge| |git3moji badge|

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

Key Features
=============

**Add Sequences from Files or URLs:**  
Easily import sequences into a SeqBank database from local files or external URLs. SeqBank supports multiple formats, enabling you to work with a wide variety of sequence data.

**Multi-threaded Downloading:**  
When dealing with large datasets or numerous sequence URLs, SeqBank can utilize multi-threaded downloading, allowing you to download and process multiple sequences in parallel, significantly speeding up the workflow.

**Database Management:**  
SeqBank provides a range of commands to manage sequence databases:

- **Listing Sequences:** List all accessions in the SeqBank database.
- **Deleting Sequences:** Remove specific sequences from the database.
- **Counting Sequences:** Quickly determine the number of sequences stored in your SeqBank.

**Integration with RefSeq and DFam:**  
SeqBank makes it simple to download and add sequences from popular genomic databases like RefSeq and DFam. This can be done directly from the command line, allowing for seamless integration into your bioinformatics pipelines.

**Data Export:**  
Sequences stored in SeqBank can be exported to widely-used formats like FASTA, making them easily shareable and compatible with other bioinformatics tools and platforms.

**Histogram Visualization:**  
SeqBank includes built-in functionality for generating histograms of sequence lengths, providing a visual summary of the data. You can save the histogram as an image file or display it interactively for immediate analysis.

**Efficient Data Filtering:**  
SeqBank allows users to filter sequences based on specific criteria before adding them to the database, ensuring only relevant data is stored.

**Customizable Workflow:**  
With support for specifying sequence formats, controlling maximum sequence additions, and utilizing temporary directories for downloads, SeqBank provides flexibility to suit various project requirements.

Installation
============

To install the latest version from the repository, you can use this command:

.. code-block:: bash

    pip install git+https://github.com/rbturnbull/seqbank.git

.. note ::

    Soon seqbank will be able to be installed using conda and PyPI


Usage
============

TODO Explain main CLI tools

e.g. 
seqbank add

seqbank count

seqbank ls

seqbank delete


.. end-quickstart


Credits
==================================

.. start-credits

* Robert Turnbull <robert.turnbull@unimelb.edu.au>
* Rafsan Al Mamun <rafsan7238@gmail.com>

.. end-credits

