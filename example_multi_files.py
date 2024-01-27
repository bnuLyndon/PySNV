# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:58:54 2022

@author: lilab
"""

import os, time

from detect_multi_samples import process_folder, parallel_process_folder
    
from pyiSNV.utils import Configuration


#parameters
kmer_length = 21
downsample = 1
snv_limit = 0.02
nb_kernel = 1

#please change to running_dir
running_dir='/home/lee/workspace/iSNV/release'

ref_file=running_dir+'/pyiSNV/GCF_009858895.2_ASM985889v3_genomic.fna'

data_dir='/media/data/ll/iSNV_test/real/1st/'
output_dir='/home/lee/workspace/iSNV/release/output/real1_3rd/'
#data_dir = '/media/data/ll/iSNV_test/ten_samples/Samples/'
#output_dir='/home/lee/workspace/iSNV/release/output/test/'

if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

config = Configuration(kmer_length, downsample, snv_limit)

T0=time.time()
if nb_kernel > 1:
    parallel_process_folder(data_dir, ref_file, output_dir, config, nb_kernel)
else:
    stat_dict = process_folder(data_dir, ref_file, output_dir, config)

T1=time.time()
print(T1-T0)
