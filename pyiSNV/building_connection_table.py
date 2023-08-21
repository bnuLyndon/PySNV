# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:58:54 2022

@author: lilab
"""

import time
import numpy as np
from hirola import HashTable

def decode(kmer, kmer_length=21, return_base=False):
    if kmer==' ':
        return -1
    kmer=int(kmer)
    decoded_kmer=np.zeros(kmer_length,dtype=np.int64)
    for i in range(kmer_length):
        decoded_kmer[i]=kmer//(4**(kmer_length-1-i))
        kmer=kmer%(4**(kmer_length-1-i))
        
    seq=''.join(list(decoded_kmer.astype(str)))

    if return_base:
        seq=seq.replace('0','A')
        seq=seq.replace('1','C')
        seq=seq.replace('2','G')
        seq=seq.replace('3','T')

    return seq

def seq2kmer(seq, kernel, tran_table, k):
    """Converts a DNA sequence split into a list of k-mers.
    The sequences in one data set do not have to share the same length.
    Returns:
         kmer_array: a numpy array of corresponding k-mer indexes.
    """

    num_seq=seq.translate(tran_table)
    num_seq=np.array(list(num_seq),dtype=np.int64)

    convolved=np.convolve(num_seq,kernel,mode='valid')

    return convolved

def dict2table(the_dict, threshold=None, key_dtype=np.int64, value_dtype=np.int64, length=0):
    if length==0:
        length=int(len(the_dict)*2)
    table = HashTable(
        length,  # <--- Maximum size for the table 
        key_dtype,  # <--- NumPy dtype 
    )

    key_array=np.array(list(the_dict.keys()))
    value_array=np.array(list(the_dict.values()))
    if threshold is not None:
        key_array=key_array[value_array>=threshold]
        value_array=value_array[value_array>=threshold]
    #value_array=np.zeros(len(the_dict),dtype=value_dtype)
    table.add(key_array)
    #for key in the_dict.keys():
    #    idx=table.add(key_dtype(key))
    #    value_array[np.array(idx,dtype=value_dtype)]=the_dict[key]

    return table, value_array


def build_tran_table():
    tran_table={}
    for i in range(65,91):
        tran_table[chr(i)]='6'
    #for i in range(97,123):
    #    tran_table[chr(i)]='4'
    tran_table['A']='0'
    tran_table['C']='1'
    tran_table['G']='2'
    tran_table['T']='3'
    return str.maketrans(tran_table)


def batch_mat_unique(mat, nb_batch):
    if nb_batch==1:
        unique_long_kmers = np.unique(mat, axis=0)
        return unique_long_kmers
    
    #t0=time.time()
    batch_size=len(mat)//nb_batch

    long_kmers_list=[]
    for i in range(nb_batch-1):
        unique_long_kmers = np.unique(mat[i*batch_size:(i+1)*batch_size], axis=0)
        long_kmers_list.append(unique_long_kmers)
        #print(i, time.time()-t0)
    unique_long_kmers = np.unique(mat[i*batch_size:(i+1)*batch_size], axis=0)
    long_kmers_list.append(unique_long_kmers)
    
    return np.vstack(long_kmers_list)
    

def build_connection_table(high_frequency_kmer_dict, ref_kmer_dict_f, kmers, \
                        kmer_length, snp_limit=0.04):

    t0=time.time()
   
    
    hfk_table, _ = dict2table(high_frequency_kmer_dict)
    rk_table, _ = dict2table(ref_kmer_dict_f)
    
    print('kmer table built', time.time()-t0)
    
    hf_mask=hfk_table.contains(kmers)
    print('high freq kmer mask built', time.time()-t0)
    r_mask=rk_table.contains(kmers)
    print('read kmer mask built', time.time()-t0)
    
    '''
    start_mask=r_mask[:-kmer_length] & ~r_mask[kmer_length:]
    end_mask=~r_mask[:-kmer_length] & r_mask[kmer_length:]
    mid_mask=~r_mask[:-kmer_length] & ~r_mask[kmer_length:]
    
    index=np.arange(len(start_mask))[start_mask]
    unique_kmers, indexes, counts = np.unique(kmers[index+kmer_length], True, True)
    '''
    
    del hfk_table, rk_table

    
    '''
    start_mask=r_mask[:-1] & ~r_mask[1:]
    start_kmers=kmers[1:][start_mask]
    unique_start_kmers, start_counts = np.unique(start_kmers, False, True)

    #long_start_mask=start_mask[:1-kmer_length] & ~r_mask[kmer_length:]
    long_start_kmer_mat=np.zeros([np.sum(start_mask),2],dtype=np.int64)
    long_start_kmer_mat[:,0]=kmers[1:-kmer_length][start_mask[:-kmer_length]]
    long_start_kmer_mat[:,1]=kmers[kmer_length:-1][start_mask[:-kmer_length]]
    unique_long_start_kmers = np.unique(long_start_kmer_mat, axis=0)
    
    
    end_mask=~r_mask[:-1] & r_mask[1:]
    end_kmers=kmers[:-1][end_mask]
    unique_end_kmers, end_counts = np.unique(end_kmers, False, True)

    long_end_kmer_mat=np.zeros([np.sum(end_mask),2],dtype=np.int64)
    long_end_kmer_mat[:,0]=kmers[1:-kmer_length][end_mask[kmer_length:]]
    long_end_kmer_mat[:,1]=kmers[kmer_length:-1][end_mask[kmer_length:]]
    unique_long_end_kmers = np.unique(long_end_kmer_mat, axis=0)
    '''
    
    hf_r_mask = hf_mask | r_mask
    connection_mask = (hf_r_mask[:-1] & hf_mask[1:]) | (hf_mask[:-1] & hf_r_mask[1:])
    #connection_mask = hf_r_mask[:-2] & hf_mask[1:-1] & hf_r_mask[2:]
    
    long_kmer_mask = connection_mask[:1-kmer_length] & hf_r_mask[kmer_length:]

    
    del hf_mask, r_mask, hf_r_mask
    
    
    index=np.arange(len(long_kmer_mask))[long_kmer_mask]
    long_kmer_mat=np.zeros([len(index),2],dtype=np.int64)
    long_kmer_mat[:,0]=kmers[index]
    long_kmer_mat[:,1]=kmers[index+kmer_length]
    
    #unique_long_kmers = np.unique(long_kmer_mat[::downsample], axis=0)
    unique_long_kmers = batch_mat_unique(long_kmer_mat, 10)
    #unique_long_kmers = np.unique(unique_long_kmers, axis=0)
    '''
    for i in range(len(unique_long_kmers)):
        loc1=ref_kmer_dict_f.get(unique_long_kmers[i,0])
        loc2=ref_kmer_dict_f.get(unique_long_kmers[i,1])
        if loc1: 
            if counts[i] < reads_abun[loc1,0]*snp_limit*0.9:
                counts[i]=0
        elif loc2:
            if counts[i] < reads_abun[loc2,0]*snp_limit*0.9:
                counts[i]=0
        else:
            continue
            if counts[i] < median_abun_threshold:
                counts[i]=0
            
    unique_long_kmers=unique_long_kmers[counts!=0]
    '''
    
    #unique_long_kmers=unique_long_kmers[counts>=max(10, lower_bound_threshold)]

    del kmers, connection_mask, long_kmer_mask, index
    
    connection_array=np.left_shift(unique_long_kmers[:,0],2)+\
        np.right_shift(unique_long_kmers[:,1],2*(kmer_length-1))
    

    unique_connections = np.unique(connection_array)

    connection_mat=np.zeros([len(unique_connections),2],dtype=np.int64)
    connection_mat[:,0]=np.right_shift(unique_connections,2)
    connection_mat[:,1]=unique_connections-\
        np.left_shift(np.right_shift(unique_connections,2*kmer_length), 2*kmer_length)

    
    print('unique connection selected', time.time()-t0)
    print(len(unique_connections), len(unique_long_kmers))


    return connection_mat, unique_long_kmers
    

if __name__ == '__main__':
    import os, psutil, time

    default_kmer_length=21
    default_chunksize = 1000000

    high_frequency_kmer_array=np.load('/home/liliandong/workspace/iSNV/temp/high_frequency_kmer_array.npy')
    exp_high_frequency_kmer_dict=dict(high_frequency_kmer_array)

    ref_db_array=np.load('/home/liliandong/workspace/iSNV/ref_db_array_f.npy')
    exp_ref_kmer_dict_f = dict(ref_db_array)

    R1_sequence_file='/home/liliandong/workspace/iSNV/test/SampleP_150_R1.fa'
    R2_sequence_file='/home/liliandong/workspace/iSNV/test/SampleP_150_R2.fa'

    T0=time.time()
    connection_table = build_connection_table(exp_high_frequency_kmer_dict, exp_ref_kmer_dict_f, \
        R1_sequence_file, R2_sequence_file, default_chunksize, default_kmer_length)
    T1=time.time()

    #open('temp/connection_dict.pickle','wb').write(pickle.dumps(connection_dict))
    
    print('time using: ', T1-T0)
    print(u'RAM using %.4f GB' % (psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024) )