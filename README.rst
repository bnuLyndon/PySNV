iSNV Detection
==============

PySNV is a tool for detecting iSNVs (intra-host variations) in tNGS sequencing data.

Install
------------

.. code-block:: sh

    git clone https://github.com/bnuLyndon/PySNV.git
    cd PySNV
    pip install .

Check Installation
------------

.. code-block:: sh

    pyisnv --help

This will display the usage information:

.. code-block:: text

    Usage: pyisnv [OPTIONS] COMMAND [ARGS]...
    PySNV is a tool for detecting iSNVs (intra-host variations) in tNGS
    sequencing data.

    Options:
    --version  Show the version and exit.
    --help     Show this message and exit.

    Commands:
    detect-multi-samples  Detect iSNVs of multiple samples using PySNV.
    detect-sample         Detect iSNVs using PySNV.

To check the installed version, run:

.. code-block:: sh

    pyisnv --version

This should output the current version:

.. code-block:: text

    pyisnv, version 1.0.0


Usage
-----

Example Files

    example.py: Demonstrates basic usage of the PySNV tool.

    example_multi_files.py: Illustrates how to process multiple input files simultaneously. Paired samples hould include '_R1' and '_R2' in the file names. Adjust the number of kernels to use parallel processing.

Note: Before running the examples, make sure to set the necessary parameters in the file.

.. code-block:: sh

    python example.py <path-to-pyiSNV> <path-to-R1-fastq>
    python example_multi_files.py <path-to-pyiSNV> <path-to-fastqs-folder>

Bash Usage
------------------

Command

.. code-block:: sh

    pyisnv detect-sample \
        --sample1 <path-to-sample1.fastq> \
        --sample2 <path-to-sample2.fastq> \
        --reference <path-to-reference_genome.fa> \
        --output <output_filename>

    pyisnv detect-multi-samples \
        --folder <path-to-samples-folder>  \
        --reference <path-to-reference_genome.fa> \
        --output <output-folder>

Parameter

``--sample1``: Path to single-end sample or first paired-end sample. Fasta, fastq and gz files are supported.

``--sample2`` (optional): Path to the second paired-end sample.

``--reference`` (required): Path to the reference genome.

``--output``: Output file path and name (default is current working directory, using sample name as output filename).

Additional Parameters

    ``--threshold``: Detection Threshold (Default: 0.02)
        The recommended detection threshold should be lager than sequencing error rate.\
    ``--kmer_length``: Kmer Length (Default: 21)
        This parameter specifies the length of k-mers to be considered during the analysis and must be an odd number smaller than 31. Kmer length should be set to ensure no duplicate kmer exists on the genome.\
    ``--downsample``: Downsample Factor (Default: 1)
        Downsamping of sequencing reads could enhance detection speed and reduce RAM usage, which could be used for high-depth sequencing samples. \
    ``--error_rate``: Sequencing Error Rate (Default: 0.01)
        Used to filter out possible false positive detection.\
    ``--indel_limit``: Maximum Indel Length (Default: 300)
        To mitigate false positive indels, especially in the case of challenging long insertions and potential impacts on estimated sequencing depths due to long deletions, a default maximum indel length of 300 is set. The recommended length threshold is 2*average_read_length.\

Cite
-----
If these codes are helpful for you, please cite our paper:

Li, L., Fu, H., Ma, W., & Li, M. (2024). PySNV for complex intra-host variation detection. Bioinformatics, btae116.

https://doi.org/10.1093/bioinformatics/btae116

https://academic.oup.com/bioinformatics/article/40/3/btae116/7616992
