# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:58:54 2022

@author: lilab
"""

import os, time
import argparse
from multiprocessing import Pool
    
from pyiSNV.utils import Configuration, get_file_name, get_file_type, accept_type_list
from detect_sample import detect_variant

def is_second(file):
    if '_R2' in file:
        return True
    
def find_paired(file):
    if '_R1' in file:
        return file.replace('_R1', '_R2')
    #elif '_1.' in file:
    #    return file.replace('_1.', '_2.')
    return ''

def process_folder(data_dir, ref_file, output_dir, config):
    files=os.listdir(data_dir)
    stat_dict = dict()
    for file in files:
        if is_second(file):
            continue
        if get_file_type(file) not in accept_type_list:
            continue
        R1_file = data_dir+file
        R2_file = find_paired(R1_file) 
        output_file = output_dir + get_file_name(R1_file) +'.csv'
        
        reads_kmer_count, genome_kmer_count, genome_kmer_rate, excluded_kmer_count,\
            excluded_kmer_rate, hf_kmer_rate, variant_kmer_rate, connection_rate,\
            utilization_rate = detect_variant(R1_file, R2_file, ref_file, config, output_file)
            
        #stat_dict[file] = [reads_kmer_count, genome_kmer_count, genome_kmer_rate, \
        #    excluded_kmer_count, excluded_kmer_rate, hf_kmer_rate, variant_kmer_rate, \
        #    connection_rate,utilization_rate]
            
    return stat_dict

def parallel_process_folder(data_dir, ref_file, output_dir, config, 
                            max_process, interval=1):

    pool = Pool(max_process)
        
    files=os.listdir(data_dir)
    for file in files:
        if is_second(file):
            continue
        if get_file_type(file) not in accept_type_list:
            continue
        R1_file = data_dir+file
        R2_file = find_paired(R1_file) 
        output_file = output_dir + get_file_name(R1_file) +'.csv'
        pool.apply_async(func=detect_variant, args=(R1_file, R2_file, ref_file, config, output_file))
        time.sleep(interval)
        #try:
        #    process_file(filepath,running_dir,temp_dir)
        #except:
        #    pass
            
    pool.close()
    pool.join()
    return True
    

    
def main():
    parser = argparse.ArgumentParser(description='Detect iSNVs using PySNV.')
    parser.add_argument('--folder', type=str, help='Path to sample folder')
    parser.add_argument('--reference', type=str, required=True, help='Path to reference genome')
    parser.add_argument('--output', type=str, default='', help='Output folder')
    
    parser.add_argument('--kernel', type=int, default=1, help='Number of kernels')
    parser.add_argument('--threshold', type=float, default=0.02, help='Detection threshold')
    parser.add_argument('--kmer_length', type=int, default=21, help='Kmer length')
    parser.add_argument('--downsample', type=int, default=1, help='Downsample factor')
    parser.add_argument('--error_rate', type=float, default=0.01, help='Sequencing error rate')
    parser.add_argument('--indel_limit', type=int, default=300, help='Maximum indel length')
    args = parser.parse_args()

    try:
        config = Configuration(kmer_length=args.kmer_length, downsample=args.downsample, 
                               snv_limit=args.threshold, error_rate=args.error_rate,
                               indel_limit=args.indel_limit)
        
        if args.kernel > 1:
            parallel_process_folder(args.folder, args.reference, args.output, config, args.kernel)
        else:
            process_folder(args.folder, args.reference, args.output, config)
        
    except Exception as e:
        print(f"Error: {e}")
    
if __name__ == '__main__':
    
    main()
