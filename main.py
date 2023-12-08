# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:58:54 2022

@author: lilab
"""

import sys, os, time
import numpy as np

from multiprocessing import Pool

from pyiSNV import build_ref_db, build_reads_db, build_connection_table, \
    output_iSNV_area, output_isnv

def format_path(path):
    return path.replace('\r','/r').replace('\\','/')

def get_file_extension(file_path):
    _, file_extension = os.path.splitext(file_path)
    return file_extension.lower()


def process_file(R1_sequence_file,running_dir,temp_path, output_dir):
    R2_sequence_file=''
    if 'R1' in R1_sequence_file:
        R2_sequence_file=R1_sequence_file.replace('R1', 'R2')

    #default_average_kmer_count=4350.3*0.818*5
    #default_average_kmer_count=12716.26*0.818/2

    running_dir=format_path(running_dir)
    R1_sequence_file=format_path(R1_sequence_file)
    R2_sequence_file=format_path(R2_sequence_file)

    example_ref_file=running_dir+'/pyiSNV/GCF_009858895.2_ASM985889v3_genomic.fna'

    output_name=os.path.split(R1_sequence_file)[-1]
    output_name=output_name.split('.')[0]
    output_path=output_dir+output_name+'.txt'
    
    if os.path.isfile(output_path):
        print(output_name, 'already exist')
        return True

    if os.path.isfile(R2_sequence_file):
        print('pair-ended data: ', output_name)
    else:
        print('single-end data: ', output_name)
    
    exp_qc_dir=temp_path+output_name+'_review/'
    if not os.path.isdir(exp_qc_dir):
        
        os.mkdir(exp_qc_dir)

    T0 = time.time()


    ref_db_array_f, ref_db_array_r, seq = build_ref_db(example_ref_file, default_kmer_length)
    ref_db_array=np.concatenate([ref_db_array_f, ref_db_array_r])

    print('1/5 ref db building complete')
    T1=time.time()
    print('time using: ', T1-T0)
    print('-------------------------------')

    exp_ref_kmer_dict = dict(ref_db_array)
    exp_ref_kmer_f_dict = dict(ref_db_array_f)
    exp_ref_kmer_r_dict = dict(ref_db_array_r)

    exp_high_frequency_kmer_dict, reads_kmers, reads_abun, com_abun, snp_list, border_list = \
        build_reads_db(R1_sequence_file, R2_sequence_file, seq, \
        exp_ref_kmer_dict, exp_ref_kmer_f_dict, exp_ref_kmer_r_dict, \
        default_kmer_length, default_snp_limit, default_downsample)
        
    print('2/5 reads db building complete')
    T1=time.time()
    print('time using: ', T1-T0)
    print('-------------------------------')

        
    connection_mat, unique_long_kmers, long_kmer_counts \
        = build_connection_table(exp_high_frequency_kmer_dict, exp_ref_kmer_f_dict, \
            reads_kmers, default_kmer_length, default_snp_limit)
    #print(len(exp_connection_dict))


    print('3/5 connection table building complete')
    T1 = time.time()
    print('time using: ', T1-T0)
    print('-------------------------------')

    
    iSNV_area_list = output_iSNV_area(seq, exp_high_frequency_kmer_dict, exp_ref_kmer_dict, \
        exp_ref_kmer_f_dict, connection_mat, unique_long_kmers, long_kmer_counts,\
        reads_abun+com_abun, default_kmer_length, default_snp_limit)
    
    print('4/5 iSNV position calculated')
    T1 = time.time()
    print('time using: ', T1-T0)
    print('-------------------------------')
    
    
    nb_iSNV, index_dict, iSNV_dict = output_isnv(seq, iSNV_area_list, snp_list, border_list, exp_qc_dir, exp_high_frequency_kmer_dict, \
        output_path, default_kmer_length, reads_abun,\
        default_snp_limit, default_downsample, default_count_limit, \
        default_average_reads_length, default_error_rate, default_vcf_threshold)
    
    print('5/5 iSNV table printed')
    print('-------------------------------')
    
    print('total types of iSNV:', nb_iSNV)
    T1=time.time()
    print('time using: ', T1-T0)
    
    return True

def process_folder(data_dir, running_dir, temp_dir, output_dir):
    files = os.listdir(data_dir)
    for file in files:
        if '_R2' in file:
            continue
        filepath=data_dir+file
        extension = get_file_extension(filepath)
        if extension in file_formats:
            process_file(filepath, running_dir, temp_dir, output_dir)
        #try:
        #    process_file(filepath, running_dir, temp_dir)
        #except:
        #    pass
            
    return True

def parallel_process_folder(data_dir,running_dir,temp_dir, output_dir, max_process=5, interval=5):

    pool = Pool(max_process)
        
    files=os.listdir(data_dir)
    for file in files:
        if '_R2' in file:
            continue
        filepath=data_dir+file
        extension = get_file_extension(filepath)
        if extension in file_formats:
            pool.apply_async(func=process_file,args=(filepath, running_dir, temp_dir, output_dir))
        #time.sleep(interval)
        #try:
        #    process_file(filepath, running_dir, temp_dir)
        #except:
        #    pass
            
    pool.close()
    pool.join()
    return True
    

#parameters
default_kmer_length=21
default_downsample=1
default_average_reads_length=150

#usually no need to change
default_snp_limit=0.02
default_count_limit=10
default_error_rate=0.01
default_vcf_threshold=0.9999

file_formats = set(['.fa', '.fq'])

#please change to running_dir
running_dir = '.'

nb_kernels=1

data_dir = 'tests/'
temp_dir = 'temp/'
output_dir = 'output/'
if nb_kernels==1:
    process_folder(data_dir, running_dir, temp_dir, output_dir) #single kernel
else:
    parallel_process_folder(data_dir, running_dir, temp_dir, output_dir) #multi kernel

