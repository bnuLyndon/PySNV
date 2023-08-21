# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:58:54 2022

@author: lilab
"""

import sys, os, time
import gzip
import numpy as np
import scipy.signal as signal

from multiprocessing import Process, Pool

from pyiSNV import build_ref_db, build_reads_db, build_connection_table, \
    output_iSNV_area, recognise_isnv

def format_path(path):
    return path.replace('\r','/r').replace('\\','/')


def process_file(R1_sequence_file,running_dir,temp_path, output_dir):
    R2_sequence_file=''

    #default_average_kmer_count=4350.3*0.818*5
    #default_average_kmer_count=12716.26*0.818/2

    running_dir=format_path(running_dir)
    R1_sequence_file=format_path(R1_sequence_file)
    R2_sequence_file=format_path(R2_sequence_file)

    example_ref_file=running_dir+'/pyiSNV/GCF_009858895.2_ASM985889v3_genomic.fna'

    output_name=os.path.split(R1_sequence_file)[-1]
    output_name=output_name.split('.')[0]
    output_path=output_dir+output_name+'.txt'

    if R2_sequence_file:
        print('pair-ended data: ', output_name)
    else:
        print('single-end data: ', output_name)

    exp_iSNV_dir=temp_path+output_name+'_kmers/'
    if os.path.isdir(exp_iSNV_dir):
        return True
        os.removedirs(exp_iSNV_dir)
    os.mkdir(exp_iSNV_dir)

    exp_qc_iSNV_dir=temp_path+output_name+'_qc/'
    if os.path.isdir(exp_qc_iSNV_dir):
        return True
        os.removedirs(exp_qc_iSNV_dir)
    os.mkdir(exp_qc_iSNV_dir)

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

    high_frequency_kmer_array, reads_ref_kmer_array, reads_kmers = build_reads_db(exp_ref_kmer_dict, \
        R1_sequence_file, R2_sequence_file, default_chunksize, \
        default_kmer_length, default_snp_limit, default_kmer_count_threshold)
        
    print('2/5 reads db building complete')
    T1=time.time()
    print('time using: ', T1-T0)
    print('-------------------------------')

        
    reads_abun=np.zeros([30000,2])
    
    for i in range(len(reads_ref_kmer_array)):
        idx=exp_ref_kmer_f_dict.get(reads_ref_kmer_array[i,0],-1)
        if idx:
            reads_abun[idx,0]=reads_ref_kmer_array[i,1]
        idx=exp_ref_kmer_r_dict.get(reads_ref_kmer_array[i,0],-1)
        if idx:
            reads_abun[idx,1]=reads_ref_kmer_array[i,1]

    average_kmer_count=np.median(reads_ref_kmer_array[:,1])
    lower_bound_threshold=default_snp_limit*0.5*average_kmer_count
    
    exp_high_frequency_kmer_dict=dict(high_frequency_kmer_array)
    
    exp_connection_dict, long_kmer_set = build_connection_table(exp_high_frequency_kmer_dict, exp_ref_kmer_f_dict, \
        reads_kmers, default_chunksize, default_kmer_length, lower_bound_threshold)
    #print(len(exp_connection_dict))


    print('3/5 connection table building complete')
    T1 = time.time()
    print('time using: ', T1-T0)
    print('-------------------------------')

    
    output_iSNV_area(seq, exp_high_frequency_kmer_dict, exp_ref_kmer_dict, \
        exp_ref_kmer_f_dict, exp_connection_dict, long_kmer_set, temp_path, output_name, default_kmer_length)
    
    print('4/5 iSNV position calculated')
    T1 = time.time()
    print('time using: ', T1-T0)
    print('-------------------------------')
    
    
    nb_iSNV = recognise_isnv(exp_iSNV_dir, exp_qc_iSNV_dir, exp_high_frequency_kmer_dict, \
        output_path, default_kmer_length, reads_abun)
    
    print('5/5 iSNV table printed')
    print('-------------------------------')
    
    print('total types of iSNV:', nb_iSNV)
    T1=time.time()
    print('time using: ', T1-T0)
    
    return True

def process_folder(data_dir, running_dir, temp_dir, output_dir):
    files = os.listdir(data_dir)
    for file in files:
        filepath=data_dir+file
        process_file(filepath, running_dir, temp_dir, output_dir)
        #try:
        #    process_file(filepath, running_dir, temp_dir)
        #except:
        #    pass
            
    return True

def multi_process_folder(data_dir,running_dir,temp_dir, output_dir, max_process, interval=5):

    pool = Pool(max_process)
        
    files=os.listdir(data_dir)
    for file in files:
        
        filepath=data_dir+file
        pool.apply_async(func=process_file,args=(filepath, running_dir, temp_dir, output_dir))
        #time.sleep(interval)
        #try:
        #    process_file(filepath, running_dir, temp_dir)
        #except:
        #    pass
            
    pool.close()
    pool.join()
    return True
    

default_kmer_length = 21
default_chunksize = 5000000

default_snp_limit = 0.02
default_kmer_count_threshold = -1

#please change to running_dir
running_dir = '.'

data_dir = 'test'
temp_dir = 'temp/'
output_dir = 'output/'
#process_folder(data_dir, running_dir, temp_dir, output_dir) #single kernel
#multi_process_folder(data_dir, running_dir, temp_dir, output_dir) #multi kernel

