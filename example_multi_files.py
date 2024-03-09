# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:58:54 2022

@author: lilab
"""

import os
import time
import sys
from pyiSNV.pyisnv import process_folder, parallel_process_folder
from pyiSNV.utils import Configuration


# parameters
kmer_length = 21
downsample = 1
snv_limit = 0.02
nb_kernel = 1

running_dir = sys.argv[1]
data_dir = sys.argv[2]
ref_file = os.path.join(running_dir,"data","GCF_009858895.2_ASM985889v3_genomic.fna")
output_dir = "output_example_multi_files"

if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

config = Configuration(kmer_length, downsample, snv_limit)

T0 = time.time()
if nb_kernel > 1:
    parallel_process_folder(data_dir, ref_file, output_dir, config, nb_kernel)
else:
    stat_dict = process_folder(data_dir, ref_file, output_dir, config)

T1 = time.time()
print(T1 - T0)
