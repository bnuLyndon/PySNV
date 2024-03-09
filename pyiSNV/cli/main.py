# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:58:54 2022

@author: lilab
"""

import sys
import click
from pyiSNV import __version__
from pyiSNV.pyisnv import detect_variant, parallel_process_folder, process_folder
from pyiSNV.utils import Configuration


@click.group()
@click.version_option(__version__)
def cli():
    "PySNV is a tool for detecting iSNVs (intra-host variations) in tNGS sequencing data."
    pass


@cli.command()
@click.option(
    "--sample1",
    type=str,
    required=True,
    help="Path to single-end fastq file or R1 paired-end fastq file",
)
@click.option(
    "--sample2",
    type=str,
    help="Path R2 paired-end fastq file (optional)",
)
@click.option(
    "--reference",
    type=str,
    required=True,
    help="Path to reference genome",
)
@click.option(
    "--output",
    type=str,
    default="output.tsv",
    help="Output file name",
)
@click.option(
    "--threshold",
    type=float,
    default=0.2,
    help="Detection threshold",
)
@click.option(
    "--kmer_length",
    type=int,
    default=21,
    help="Kmer length",
)
@click.option(
    "--downsample",
    type=int,
    default=1,
    help="Downsample factor",
)
@click.option(
    "--error_rate",
    type=float,
    default=0.1,
    help="Sequencing error rate",
)
@click.option(
    "--indel_limit",
    type=int,
    default=300,
    help="Maximum indel length",
)
def detect_sample(
    sample1,
    sample2,
    reference,
    output,
    threshold,
    kmer_length,
    downsample,
    error_rate,
    indel_limit,
):
    """
    Detect iSNVs using PySNV.
    """
    try:
        config = Configuration(
            kmer_length=kmer_length,
            downsample=downsample,
            snv_limit=threshold,
            error_rate=error_rate,
            indel_limit=indel_limit,
        )
        detect_variant(sample1, sample2, reference, config, output)
    except Exception as err:
        click.secho(err, err=True, fg="red")
        sys.exit(1)


@cli.command()
@click.option(
    "--folder",
    type=str,
    required=True,
    help="Path to sample folder",
)
@click.option(
    "--reference",
    type=str,
    required=True,
    help="Path to reference genome",
)
@click.option(
    "--output",
    type=str,
    default="output",
    help="Output folder",
)
@click.option(
    "--kernel",
    type=int,
    default=1,
    help="Number of kernels",
)
@click.option(
    "--threshold",
    type=float,
    default=0.2,
    help="Detection threshold",
)
@click.option(
    "--kmer_length",
    type=int,
    default=21,
    help="Kmer length",
)
@click.option(
    "--downsample",
    type=int,
    default=1,
    help="Downsample factor",
)
@click.option(
    "--error_rate",
    type=float,
    default=0.1,
    help="Sequencing error rate",
)
@click.option(
    "--indel_limit",
    type=int,
    default=300,
    help="Maximum indel length",
)
def detect_multi_samples(
    folder,
    reference,
    output,
    kernel,
    threshold,
    kmer_length,
    downsample,
    error_rate,
    indel_limit,
):
    """
    Detect iSNVs of multiple samples using PySNV.
    """
    try:
        config = Configuration(
            kmer_length=kmer_length,
            downsample=downsample,
            snv_limit=threshold,
            error_rate=error_rate,
            indel_limit=indel_limit,
        )

        if kernel > 1:
            parallel_process_folder(folder, reference, output, config, kernel)
        else:
            process_folder(folder, reference, output, config)

    except Exception as err:
        click.secho(err, err=True, fg="red")
        sys.exit(1)
