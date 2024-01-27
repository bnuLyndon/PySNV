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

Usage
-----

1. Modify the `running_dir` variable in `main.py` to set the running directory.
2. Set the sample file(s) and genome file in `main.py`.
3. Adjust the number of kernels in `main.py`.
4. Run `main.py`.

Main.py Parameters
------------------

```bash
python main.py --sample1 path/to/sample1.fastq --sample2 path/to/sample2.fastq --reference path/to/reference_genome.fa --output output_filename

    --sample1: Path to single-end sample or first paired-end sample.
    --sample2 (optional): Path to the second paired-end sample.
    --reference (required): Path to the reference genome.
    --output: Output file name (default is current working directory).

Additional Parameters

    --threshold: Detection threshold (default: 0.02).
    --kmer_length: Kmer length (default: 21).
    --downsample: Downsample factor (default: 1).
    --error_rate: Sequencing error rate (default: 0.01).
    --indel_limit: Maximum indel length (default: 300).

Example Files

    example.py: Demonstrates basic usage of the PySNV tool.
    detect_sample.py: Example script showcasing the detection process for a single sample.
    example_multi_files.py: Illustrates how to process multiple input files simultaneously.

Note: Before running the examples, make sure to set the necessary parameters in the main.py file.

python

# example.py
# ... (example script content)

# detect_sample.py
# ... (example script content)

# example_multi_files.py
# ... (example script content)
