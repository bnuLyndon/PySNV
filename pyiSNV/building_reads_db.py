# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:58:54 2022

@author: lilab
"""

import os
import psutil
import time
import pandas as pd
import numpy as np
import gzip

from hirola import HashTable

from pyiSNV.utils import decode, tran_table, dict2table

class SampleMapBuilder:
    def __init__(self, ref_kmer_dict, ref_kmer_f_dict, ref_kmer_r_dict, \
                       kmer_length = 21, snv_limit = 0.04, downsample = 2):        
        self.ref_kmer_dict = ref_kmer_dict
        self.ref_kmer_f_dict = ref_kmer_f_dict
        self.ref_kmer_r_dict = ref_kmer_r_dict
        
        self.kmer_length = kmer_length
        self.snv_limit = snv_limit
        self.downsample = downsample
        
    def build_map(self, R1_file, R2_file = ''):
        return build_reads_db(R1_file, R2_file, self.seq, self.ref_kmer_dict,  
                              self.ref_kmer_f_dict, self.ref_kmer_r_dict, \
                           self.kmer_length, self.snv_limit, self.downsample)

def reads2kmer(ref_kmer_dict, ref_kmer_f_dict, ref_kmer_r_dict, R1_file, \
               kmer_length, downsample):
    t0 = time.time()
    kernel = 4**np.array(range(kmer_length),dtype = np.int64)

    '''
    ref_file = '/home/liliandong/workspace/iSNV/DB/GCF_009858895.2_ASM985889v3_genomic.fna'
    with open(ref_file, 'r') as handle:
        for rec in SeqIO.parse(handle, 'fasta'):
            seq = ''.join(rec.seq)
            kmers = seq2kmer(seq[:-33], kernel, tran_table, kmer_length)

    print(len(kmers),len(np.unique(kmers)))
    unique_kmers, index, counts = np.unique(kmers,return_index = True,return_counts = True)
    '''

    shift_size = kmer_length//2
    
    case_name = os.path.split(R1_file)[-1]
    case_name = case_name.split('.')[0]
    
    
    if R1_file[-3:] == '.gz':
        f = gzip.GzipFile(R1_file)
        lines = f.readlines()
    else:
        with open(R1_file,'r') as f:
            lines = f.readlines()
            
    print(case_name, 'loaded, using', time.time()-t0)
            
    if lines[0][0] == '>' or lines[0][0] == 62:
        lines = lines[1::2]
    else:
        lines = lines[1::4]

    reads_string = ''#.join(lines)

    df = pd.DataFrame(lines[::downsample])
    del lines
    
    if R1_file[-3:] == '.gz':
        df = df[0].str.decode('utf-8')
        numseq_list = df.str.translate(tran_table)
    else:
        numseq_list = df[0].str.translate(tran_table)

    del df
    num_seq = ''.join(numseq_list)
    del numseq_list
    num_seq = bytes(num_seq, encoding = 'utf8')
    num_array = np.frombuffer(num_seq, dtype = np.int8)-48
    num_array = num_array.astype(np.int16)
    del num_seq
    convolved = np.convolve(num_array,kernel,mode = 'same')
    
    rev_num_seq = np.subtract(3,num_array)
    rev_num_seq = np.abs(rev_num_seq)
    rev_convolved = np.convolve(rev_num_seq,kernel[::-1],mode = 'same')
    
    fore_after_ind = np.where(num_array[:-shift_size] == 4)[0]
    del num_array, rev_num_seq
        #temp_kmers = convolved[minimizer_count[i]:minimizer_count[i+1]]
    for i in range(kmer_length):
        convolved[fore_after_ind[:-1]-i+shift_size] = -1
        rev_convolved[fore_after_ind[:-1]-i+shift_size] = -1
        
    convolved[-shift_size:] = -1
    convolved[:shift_size] = -1
    rev_convolved[-shift_size:] = -1
    rev_convolved[:shift_size] = -1
    
    print('kmer extracted, using', time.time()-t0)
    
    frk_table, _ = dict2table(ref_kmer_f_dict, value_dtype = np.int32)
    
    r_mask = frk_table.contains(convolved)
    filtered_convolved = np.zeros(convolved.shape,int)
    ind = 0
    for i in range(len(fore_after_ind)-1):
        if True in r_mask[fore_after_ind[i]:fore_after_ind[i+1]]:
            length = fore_after_ind[i+1]-fore_after_ind[i]
            filtered_convolved[ind:ind+length] = convolved[fore_after_ind[i]:fore_after_ind[i+1]]
            ind += length
    filtered_convolved = filtered_convolved[:ind]
    del convolved
    
    print('forward reads kmer extracted, using', time.time()-t0)
            
    rev_convolved = rev_convolved[::-1]
    r_mask = frk_table.contains(rev_convolved)
    filtered_rev_convolved = np.zeros(rev_convolved.shape,int)
    ind = 0
    for i in range(len(fore_after_ind)-1):
        if True in r_mask[fore_after_ind[i]:fore_after_ind[i+1]]:
            length = fore_after_ind[i+1]-fore_after_ind[i]
            filtered_rev_convolved[ind:ind+length] = rev_convolved[fore_after_ind[i]:fore_after_ind[i+1]]
            ind += length
    filtered_rev_convolved = filtered_rev_convolved[:ind]
    del rev_convolved
    
    convolved = np.hstack([filtered_convolved, filtered_rev_convolved])
    del filtered_convolved, filtered_rev_convolved
    
    print('reverse reads kmer extracted, using', time.time()-t0)

    rk_table, _ = dict2table(ref_kmer_dict, value_dtype = np.int32)
    r_mask = rk_table.contains(convolved)
    #rk_pos = rk_index_array[r_index]
    #rk_pos[r_index =  = -1] = -1
    

    print('kmers split, using', time.time()-t0)
    
    
    unique_kmers, counts = np.unique(convolved[r_mask], return_counts = True)
    reads_abun_dict = dict(zip(unique_kmers, counts))
    
    reads_abun = np.zeros(30000,dtype = int)
    for key in reads_abun_dict.keys():
        reads_abun[ref_kmer_dict.get(key)] += reads_abun_dict[key]
    
    print('abundance estimated, using', time.time()-t0)
    
    print(u'RAM using %.4f GB' % (psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024) )
    
    return convolved, r_mask, reads_abun, reads_string


def get_masked_kmer(kmers, shift_size, window_size = 1):
    
    kmers_high_bits = (kmers>>2*(shift_size+window_size))<<2*(shift_size+window_size)
    kmers_low_bits = kmers-((kmers>>2*(shift_size))<<2*(shift_size))

    return kmers_high_bits+kmers_low_bits


def get_snp_kmer(hf_kmers, hf_kmers_counts, ref_kmer_f_dict, reads_abun, window_size, kmer_length, snv_limit):
    
    snp_list = []
    hf_abun = np.zeros(reads_abun.shape,int)
    
    shift_size = kmer_length//2
    
    ref_kmers = np.array(list(ref_kmer_f_dict.keys()))
    ref_kmer_loc = np.array(list(ref_kmer_f_dict.values()))
    
    masked_ref_kmers = get_masked_kmer(ref_kmers, shift_size)
    masked_ref_kmer_dict = dict(zip(masked_ref_kmers, ref_kmer_loc))
    
    loc2kmer_dict = dict(zip(ref_kmer_loc, ref_kmers))
    
    masked_hf_kmers = get_masked_kmer(hf_kmers, shift_size)
    
    for i in range(len(masked_hf_kmers)):
        masked_kmer = masked_hf_kmers[i]
        if masked_ref_kmer_dict.get(masked_kmer):
            loc = masked_ref_kmer_dict.get(masked_kmer)
            #hf_abun[loc-shift_size:loc+shift_size+1]+ = hf_kmers_counts[i]
            hf_abun[loc] += hf_kmers_counts[i]
            if hf_kmers_counts[i]>snv_limit*(reads_abun[loc]+hf_kmers_counts[i]):
                snp_list.append([loc, decode(loc2kmer_dict[loc]), decode(hf_kmers[i]), \
                                 hf_kmers[i], hf_kmers_counts[i]])
            else:
                hf_kmers[i] = -1
    
    return snp_list, hf_abun

def batch_unique(kmers, nb_batch):
    
    if nb_batch == 1:
        unique_kmers, counts = np.unique(kmers, return_counts = True)
        return unique_kmers, counts
    
    t0 = time.time()
    batch_size = len(kmers)//nb_batch

    kmer_count_list = []
    for i in range(nb_batch-1):
        unique_kmers, counts = np.unique(kmers[i*batch_size:(i+1)*batch_size], return_counts = True)
        kmer_count_list.append([unique_kmers, counts])
        #print(i, time.time()-t0)
    unique_kmers, counts = np.unique(kmers[(nb_batch-1)*batch_size:], return_counts = True)
    kmer_count_list.append([unique_kmers, counts])
    #print(nb_batch-1, time.time()-t0)
    print('got unique kmers in batches')
    
    batch_length = [len(item[0]) for item in kmer_count_list]
    table_length = sum(batch_length)
    table = HashTable(
        table_length,  # <--- Maximum size for the table 
        np.int64,  # <--- NumPy dtype 
    )
    count_array = np.zeros(table_length, dtype = np.int64)
    
    for item in kmer_count_list:
        insert_index = table.add(item[0])
        count_array[insert_index] += item[1]
        #print(table.length, time.time()-t0)
        
    return table.keys[:table.length], count_array[:table.length]

def border_process(ref_kmer_f_dict, seq, convolved, reads_abun, kmer_length):

    short_dict = dict()
    for key in ref_kmer_f_dict.keys():
        kmer = key>>2*10
        if not short_dict.get(kmer):
            short_dict[kmer] = []
        short_dict[kmer].append(ref_kmer_f_dict[key])
    
    starting_dict = dict()
    ending_dict = dict()
    border_dict = dict()
    for key in ref_kmer_f_dict.keys():
        loc = ref_kmer_f_dict[key]
        if reads_abun[loc] != 0 and reads_abun[loc+1] == 0:
            ending_dict[key] = loc
            border_dict[key] = loc
        if loc >= 1 and reads_abun[loc-1] == 0 and reads_abun[loc] != 0:
            starting_dict[key] = loc
            border_dict[key] = loc

    border_table, value_array = dict2table(border_dict)
    border_mask = border_table.contains(convolved)
    border_mask[-kmer_length:] = False
    
    index = np.arange(len(border_mask))[border_mask]
    
    ending_kmer_mat = np.zeros([len(index),kmer_length],int)
    #print(ending_kmer_mat.shape, convolved.shape, border_mask.shape)
    for i in range(kmer_length):
        ending_kmer_mat[:,i] = convolved[index+i]
    ending_kmer_mat = ending_kmer_mat[ending_kmer_mat[:,1] != -1]
    ending_kmer, ending_counts = np.unique(ending_kmer_mat, return_counts = True)
    
    ending_kmer_dict = dict(zip(ending_kmer, ending_counts))
    ending_area_dict = dict()
    for i in range(len(ending_kmer_mat)):
        temp_array = ending_kmer_mat[i,:]
        if ref_kmer_f_dict.get(temp_array[1]):
            continue
        temp_count = [ending_kmer_dict.get(item) for item in temp_array]
        if min(temp_count) == 1:
            continue
        if temp_array[0] not in ending_area_dict.keys():
            ending_area_dict[temp_array[0]] = []
        ending_area_dict[temp_array[0]].append(temp_array)
        
    ending_snv_list = []
    ending_seq_set = set()
    for key in ending_area_dict.keys():
        loc = ref_kmer_f_dict.get(key)
        for temp_array in ending_area_dict[key]:
            for j in range(1, kmer_length):
                kmer = temp_array[j]
                if kmer == -1:
                    break
                short_kmer = kmer-((kmer>>2*10)<<2*10)
                if not short_dict.get(short_kmer):
                    continue
                for item in short_dict.get(short_kmer):
                    if 10 <= item-loc < 2*kmer_length:
                        ref_seq = seq[loc: item+10+1]
                        var_seq = '{}{}'.format(decode(key), decode(kmer)[kmer_length-j:])
                        if var_seq[:kmer_length+10+1] not in ending_seq_set:                 
                            #temp_seq = [decode(elem) for elem in temp_array[:j+1]]
                            temp_count = [ending_kmer_dict.get(elem) for elem in temp_array[:j+1]]
                            ending_seq_set.add(var_seq[:kmer_length+10+1])
                            ending_snv_list.append([loc, ref_seq, var_seq, temp_array[:j+1], np.median(temp_count)])
                
                
    starting_kmer_mat = np.zeros([len(index),kmer_length],int)
    for i in range(kmer_length):
        starting_kmer_mat[:,kmer_length-1-i] = convolved[index-i]
    starting_kmer_mat = starting_kmer_mat[starting_kmer_mat[:, -2] != -1]
    
    starting_kmer, starting_counts = np.unique(starting_kmer_mat, return_counts = True)
    starting_kmer_dict = dict(zip(starting_kmer, starting_counts))
    starting_area_dict = dict()
    for i in range(len(starting_kmer_mat)):
        temp_array = starting_kmer_mat[i,:]
        if ref_kmer_f_dict.get(temp_array[-2]):
            continue
        temp_count = [starting_kmer_dict.get(item) for item in temp_array]
        if min(temp_count) == 1:
            continue
        if temp_array[-1] not in starting_area_dict.keys():
            starting_area_dict[temp_array[-1]] = []
        starting_area_dict[temp_array[-1]].append(temp_array)
        
    starting_snv_list = []
    starting_seq_set = set()
    for key in starting_area_dict.keys():
        loc = ref_kmer_f_dict.get(key)
        for temp_array in starting_area_dict[key]:
            for j in range(kmer_length-1):
                kmer = temp_array[kmer_length-2-j]
                if kmer == -1:
                    break
                short_kmer = kmer>>2*11
                if not short_dict.get(short_kmer):
                    continue
                for item in short_dict.get(short_kmer):
                    #print(loc,item,short_kmer)
                    if 10 <= loc-item < 2*kmer_length:
                        ref_seq = seq[item+1: loc+kmer_length]
                        var_seq = '{}{}'.format(decode(kmer)[:j+1], decode(key))
                        if var_seq[-10-kmer_length:] not in starting_seq_set:
                            temp_count = [starting_kmer_dict.get(elem) for elem in temp_array[kmer_length-2-j:]]
                            #temp_seq = [decode(elem) for elem in temp_array[kmer_length-2-j:]]
                            starting_seq_set.add(var_seq[-10-kmer_length:])
                            starting_snv_list.append([item+1, ref_seq, var_seq, \
                                                      temp_array[kmer_length-2-j:], np.median(temp_count)])

    return starting_snv_list+ending_snv_list

def build_reads_db(R1_file, R2_file, seq, ref_kmer_dict, ref_kmer_f_dict, ref_kmer_r_dict, \
                   kmer_length = 21, snv_limit = 0.04, downsample = 2):

    t0 = time.time()
    
    convolved, reads_kmer_mask, reads_abun, reads_string = reads2kmer(ref_kmer_dict,   
                       ref_kmer_f_dict, ref_kmer_r_dict, R1_file, kmer_length, downsample)
      
    #reads_ref_kmer_dict = dict(reads_ref_kmer_array)
    reads_string_2 = ''

    if os.path.isfile(R2_file):
        
        convolved2, reads_kmer_mask2, reads_abun_2, reads_string_2 = reads2kmer(ref_kmer_dict, 
                           ref_kmer_f_dict, ref_kmer_r_dict, R2_file, kmer_length, downsample)
        
        reads_abun = reads_abun_2+reads_abun
        convolved = np.hstack([convolved, convolved2])
        reads_kmer_mask = np.hstack([reads_kmer_mask, reads_kmer_mask2])
        
        del convolved2, reads_kmer_mask2, reads_abun_2
        
    unique_kmers, counts = batch_unique(convolved[~reads_kmer_mask], 10)
    #unique_kmers, counts = np.unique(convolved[~reads_kmer_mask], return_counts = True)
    mask = (counts>2) & (unique_kmers>2**(kmer_length//2))
    unique_kmers = unique_kmers[mask]
    counts = counts[mask]
    
    high_frequency_kmer_dict = dict(zip(unique_kmers, counts))
    
    print('kmer count calculated, using', time.time()-t0)
    
    kmer_count_threshold = np.median(reads_abun)
    
    mask = counts >= kmer_count_threshold  
    
    lf_kmer = unique_kmers[~mask]
    lf_counts = counts[~mask]

    
    unique_kmers = unique_kmers[mask]
    counts = counts[mask]
    hf_snp_list, hf_abun = get_snp_kmer(unique_kmers, counts, ref_kmer_f_dict, \
                                    reads_abun, kmer_length//2+1, kmer_length, snv_limit)

    lf_snp_list, lf_abun = get_snp_kmer(lf_kmer, lf_counts, ref_kmer_f_dict, \
                                 reads_abun+hf_abun, kmer_length//2+1, kmer_length, snv_limit)

    
    print('simple snp found, using', time.time()-t0)
    
    border_snv_list = border_process(ref_kmer_f_dict, seq, convolved, reads_abun, kmer_length)
    
    return high_frequency_kmer_dict, convolved, reads_abun, hf_abun+lf_abun, \
        hf_snp_list+lf_snp_list, border_snv_list

if __name__ == '__main__':
    import os, psutil, time

    default_chunksize = 1000000
    R1_sequence_file = '/home/liliandong/workspace/iSNV/test/SampleP_150_R1.fa'
    R2_sequence_file = '/home/liliandong/workspace/iSNV/test/SampleP_150_R2.fa'

    ref_db_array = np.load('/home/liliandong/workspace/iSNV/ref_db_array.npy')
    exp_ref_kmer_dict = dict(ref_db_array)

    default_kmer_length = 21
    default_threshold = 135*5

    T0 = time.time()
    high_frequency_kmer_array, _ = build_reads_db(exp_ref_kmer_dict, R1_sequence_file, R2_sequence_file, \
        default_chunksize, default_kmer_length, default_threshold)
    T1 = time.time()

    #np.save('temp/high_frequency_kmer_array.npy', high_frequency_kmer_array)
    
    print('time using: ', T1-T0)
    print(u'RAM using %.4f GB' % (psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024) )