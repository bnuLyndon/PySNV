# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:58:54 2022

@author: lilab
"""

import sys
from os import path
from pyiSNV.pyisnv import detect_variant
from pyiSNV.utils import Configuration, format_path, get_file_name


running_dir = sys.argv[1]
R1_sequence_file = sys.argv[2]
ref_file = path.join(running_dir,"data","GCF_009858895.2_ASM985889v3_genomic.fna")

if "R1" in R1_sequence_file:
    R2_sequence_file = R1_sequence_file.replace("R1", "R2")


config = Configuration(downsample=2)

pysnv_dir = format_path(running_dir)
R1_sequence_file = format_path(R1_sequence_file)
R2_sequence_file = format_path(R2_sequence_file)

output_name = get_file_name(R1_sequence_file)
output_path = path.join(pysnv_dir, "output", output_name + ".txt")

detect_variant(R1_sequence_file, R2_sequence_file, ref_file, config, output_path)
