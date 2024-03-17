#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 22:03:37 2023

@author: lee
"""

import os
import numpy as np
from hirola import HashTable
from dataclasses import dataclass


@dataclass
class Configuration:
    kmer_length: int = 21
    downsample: int = 1
    snv_limit: float = 0.02
    error_rate: float = 0.01
    indel_limit: int = 300
    verbose: bool = False
    p_threshold: float = 0.9999
    average_reads_length: int = 150
    count_limit: int = 10
    correct_factor: float = 0.8


def format_path(path):
    return path.replace("\r", "/r").replace("\\", "/")


def get_file_name(filepath):
    file_name = os.path.split(filepath)[-1]
    file_name = file_name.split(".")[0]
    return file_name


def get_file_type(file_name):
    _, file_extension = os.path.splitext(file_name)
    return file_extension


def encode(seq, kernel, tran_table):
    num_seq = seq.translate(tran_table)
    num_seq = np.array(list(num_seq), dtype=np.int64)
    convolved = np.convolve(num_seq, kernel, mode="valid")

    assert len(convolved) == 1
    return convolved[0]


def decode(kmer, kmer_length, return_base=True):
    if kmer == " ":
        return -1
    kmer = int(kmer)
    decoded_kmer = np.zeros(kmer_length, dtype=np.int64)
    for i in range(kmer_length):
        decoded_kmer[i] = kmer // (4 ** (kmer_length - 1 - i))
        kmer = kmer % (4 ** (kmer_length - 1 - i))

    seq = "".join(list(decoded_kmer.astype(str)))

    if return_base:
        seq = seq.replace("0", "A")
        seq = seq.replace("1", "C")
        seq = seq.replace("2", "G")
        seq = seq.replace("3", "T")

    return seq


def build_tran_table_old():
    tran_table = {}
    for i in range(65, 91):
        tran_table[chr(i)] = "6"

    tran_table["A"] = "0"
    tran_table["C"] = "1"
    tran_table["G"] = "2"
    tran_table["T"] = "3"
    return str.maketrans(tran_table)


def build_tran_table():
    tran_table = {}
    for i in range(65, 91):
        tran_table[chr(i)] = "4"

    tran_table["A"] = "0"
    tran_table["C"] = "1"
    tran_table["G"] = "2"
    tran_table["T"] = "3"
    tran_table["\n"] = "4"
    return str.maketrans(tran_table)


def dict2table_old(the_dict, key_dtype=np.int64, value_dtype=np.int64, length=0):
    if length == 0:
        length = int(len(the_dict) * 2)
    table = HashTable(
        length,  # <--- Maximum size for the table
        key_dtype,  # <--- NumPy dtype
    )

    table.add(np.array(list(the_dict.keys())))
    value_array = np.array(list(the_dict.values()))

    return table, value_array


def dict2table(
    the_dict, key_dtype=np.int64, value_dtype=np.int64, length=0, threshold=0
):
    if length == 0:
        length = int(len(the_dict) * 2)
    table = HashTable(
        length,  # <--- Maximum size for the table
        key_dtype,  # <--- NumPy dtype
    )

    key_array = np.array(list(the_dict.keys()))
    value_array = np.array(list(the_dict.values()))
    if threshold > 0:
        key_array = key_array[value_array >= threshold]
        value_array = value_array[value_array >= threshold]
    table.add(key_array)

    return table, value_array


def is_second(file):
    if "_R2" in file:
        return True


def find_paired(file):
    if "_R1" in file:
        return file.replace("_R1", "_R2")

    return ""


TRAN_TABLE = build_tran_table()
ACCEPT_TYPE_LIST = [".fasta", ".fastq", ".fq", ".fq", ".gz"]
