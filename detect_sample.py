# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:58:54 2022

@author: lilab
"""

import os, time
import argparse
import numpy as np

from pyiSNV import GenomeMapBuilder, SampleMapBuilder, ConnectMapBuilder, \
    VariantLocator, VariantCaller
    
from pyiSNV.utils import Configuration, format_path, get_file_name, get_file_type, accept_type_list


def detect_variant(R1_sequence_file, R2_sequence_file, ref_file, config, output_path):
    
    assert get_file_type(R1_sequence_file) in accept_type_list
    if output_path == '':
        output_path = get_file_name(R1_sequence_file)
        
    
    T0=time.time()

    ref_builder = GenomeMapBuilder(config)

    ref_builder.build_ref_db(ref_file)


    print('1/5 ref db building complete')
    T1=time.time()
    print('time using: ', np.round(T1-T0,2), 'seconds')
    print('-------------------------------')

    sample_mapper = SampleMapBuilder(config)

    high_frequency_kmer_dict, reads_kmers, reads_depth = \
        sample_mapper.build_reads_db(R1_sequence_file, R2_sequence_file, ref_builder)

    print('2/5 reads db building complete')
    T1=time.time()
    print('time using: ', np.round(T1-T0,2), 'seconds')
    print('-------------------------------')


    connection_builder = ConnectMapBuilder(config)

    connection_mat, unique_long_kmers, long_kmer_counts = \
        connection_builder.build_connection_table(high_frequency_kmer_dict, reads_kmers, ref_builder)

    print('3/5 connection table building complete')
    T1=time.time()
    print('time using: ', np.round(T1-T0,2), 'seconds')
    print('-------------------------------')
        

    variant_locator = VariantLocator(config)

    iSNV_area_list = variant_locator.output_iSNV_area(high_frequency_kmer_dict, \
        connection_mat, unique_long_kmers, long_kmer_counts, reads_depth, ref_builder)

    print('4/5 iSNV position calculated')
    T1=time.time()
    print('time using: ', np.round(T1-T0,2), 'seconds')
    print('-------------------------------')

    #return True

    variant_caller = VariantCaller(config)

    nb_iSNV = variant_caller.output_isnv(iSNV_area_list, high_frequency_kmer_dict, \
                    reads_depth, output_path, ref_builder, sample_mapper)
        
    print('total cases of iSNV:', nb_iSNV)

    print('5/5 iSNV table output')
    T1=time.time()
    print('time using:', np.round(T1-T0,2), 'seconds')
    print('-------------------------------')

    reads_kmer_count = len(reads_kmers[reads_kmers != -1]) // 2
    genome_kmer_count = sample_mapper.nb_genome_kmer // 2

    print('total kmer:', reads_kmer_count)
    print('counts of genome kmer:', genome_kmer_count)
    print('filtered-in kmer:', variant_locator.hf_kmer_count)
    print('connected high freq kmer:', variant_locator.connected_kmer_count)
    print('variant kmer:', variant_caller.connected_kmer_counts)
    print('-------------------------------')

    genome_kmer_rate = genome_kmer_count/ reads_kmer_count
    excluded_kmer_count = (reads_kmer_count - genome_kmer_count - variant_caller.connected_kmer_counts)
    excluded_kmer_rate = excluded_kmer_count / reads_kmer_count
    hf_kmer_rate = variant_locator.hf_kmer_count / reads_kmer_count
    variant_kmer_rate = variant_caller.connected_kmer_counts / reads_kmer_count
    connection_rate = min(1, variant_locator.connected_kmer_count / variant_locator.hf_kmer_count)
    utilization_rate = min(1, variant_caller.connected_kmer_counts / variant_locator.connected_kmer_count)
        
    print('genome kmer proportion:', f"{genome_kmer_rate * 100:.2f}%")
    print('excluded kmer proportion:', f"{excluded_kmer_rate * 100:.2f}%")
    print('filtered-in kmer proportion:', f"{hf_kmer_rate * 100:.2f}%")
    print('filtered-in kmer connection rate:', f"{connection_rate * 100:.2f}%")
    print('connected kmer utlization rate:', f"{utilization_rate * 100:.2f}%")
    print('variant kmer proportion:', f"{variant_kmer_rate * 100:.2f}%")

    print('-------------------------------')
    return reads_kmer_count, genome_kmer_count, genome_kmer_rate, excluded_kmer_count,\
        excluded_kmer_rate, hf_kmer_rate, variant_kmer_rate, connection_rate, utilization_rate
    
def main():
    parser = argparse.ArgumentParser(description='Detect iSNVs using PySNV.')
    parser.add_argument('--sample1', type=str, help='Path to single-end sample or first paired-end sample')
    parser.add_argument('--sample2', type=str, help='Path to second paired-end sample (optional)')
    parser.add_argument('--reference', type=str, required=True, help='Path to reference genome')
    parser.add_argument('--output', type=str, default='', help='Output file name')
    
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
        #detect_sample(R1_sequence_file, R2_sequence_file, ref_file, config, output_path)
        detect_variant(args.sample1, args.sample2, args.reference, config, args.output)
    except Exception as e:
        print(f"Error: {e}")
    
if __name__ == '__main__':
    
    main()


