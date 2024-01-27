iSNV Detection
==============

PySNV is a tool for detecting iSNVs (intra-host variations) in tNGS sequencing data.

Requirements
------------

Ensure you have the following Python libraries installed:

1. `numpy`
2. `pandas`
3. `psutil`
4. `Bio`
5. `hirola`

Python Usage
-----
1. Modify the `running_dir` variable in `example.py` to set the running directory.

2. Set the sample file(s) and genome file in `example.py`.

3. Run `example.py`.

Example Files

    example.py: Demonstrates basic usage of the PySNV tool.

    detect_sample.py: Main function of processing for a single sample.

    example_multi_files.py: Illustrates how to process multiple input files simultaneously. Paired samples hould include '_R1' and '_R2' in the file names. Adjust the number of kernels to use parallel processing.

    detect_sample.py: Main function of processing for a sample folder.

Note: Before running the examples, make sure to set the necessary parameters in the file.

Bash Usage
------------------
Command

``python detect_sample.py --sample1 path/to/sample1.fastq --sample2 path/to/sample2.fastq --reference path/to/reference_genome.fa --output output_filename``

Parameters

``--sample1``: Path to single-end sample or first paired-end sample. Fasta, fastq and gz files are supported.

``--sample2`` (optional): Path to the second paired-end sample.

``--reference`` (required): Path to the reference genome.

``--output``: Output file path and name (default is current working directory, using sample name as output filename).

Additional Parameters

    ``--threshold``: Detection Threshold (Default: 0.02)
        The recommended detection threshold should be lager than sequencing error rate.\
    ``--kmer_length``: Kmer Length (Default: 21)
        Currently capped at 30 (must be an odd number), the kmer length is automatically assessed by PySNV to ensure the absence of duplicate kmers in the genome. A slightly larger kmer length than the threshold, where duplicate kmers are absent, is recommended for optimal performance.\
    ``--downsample``: Downsample Factor (Default: 1)
        To enhance speed and reduce RAM usage, especially in high-depth sequencing scenarios, set the downsample factor (default: 1) to a value greater than or equal to 2. It is suggested to maintain a post-downsampling depth greater than 300X to preserve detection accuracy.\
    ``--error_rate``: Sequencing Error Rate (Default: 0.01)
        Used to filter out possible false positive detection..\
    -``-indel_limit``: Maximum Indel Length (Default: 300)
        To mitigate false positive indels, especially in the case of challenging long insertions and potential impacts on estimated sequencing depths due to long deletions, a default maximum indel length of 300 is set. The recommended length threshold is 2*average_read_depth.\

