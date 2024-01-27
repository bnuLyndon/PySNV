# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:58:54 2022

@author: lilab
"""

import os

from detect_sample import detect_variant
    
from pyiSNV.utils import Configuration, format_path, get_file_name


#parameters
kmer_length = 21
downsample = 2
snv_limit = 0.02

#please change to running_dir
running_dir='/home/lee/workspace/iSNV/release'

ref_file=running_dir+'/pyiSNV/GCF_009858895.2_ASM985889v3_genomic.fna'

R1_sequence_file='/home/lee/workspace/iSNV/data/SampleP150_R1.fq.gz'
#R2_sequence_file=''
if 'R1' in R1_sequence_file:
    R2_sequence_file=R1_sequence_file.replace('R1', 'R2')


config = Configuration(kmer_length, downsample, snv_limit)

running_dir=format_path(running_dir)
R1_sequence_file=format_path(R1_sequence_file)
R2_sequence_file=format_path(R2_sequence_file)

output_name=get_file_name(R1_sequence_file)
output_path=running_dir+'/output/'+output_name+'.txt'

detect_variant(R1_sequence_file, R2_sequence_file, ref_file, config, output_path)
