# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:58:54 2022

@author: lilab
"""

import numpy as np
from Bio import SeqIO

def seq2kmer(seq, kernel, tran_table, k):
    """Converts a DNA sequence split into a list of k-mers.
    The sequences in one data set do not have to share the same length.
    Returns:
         kmer_array: a numpy array of corresponding k-mer indexes.
    """

    num_seq=seq.translate(tran_table)
    num_seq=np.array(list(num_seq),dtype=np.int64)

    convolved=np.convolve(num_seq,kernel,mode='valid')

    rev_num_seq=3-num_seq
    kernel=kernel[::-1]
    rev_convolved=np.convolve(rev_num_seq, kernel, mode='valid')
    return convolved, rev_convolved

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


def build_ref_db(ref_file, kmer_length=21):

    tran_table=build_tran_table()
    kernel=4**np.array(range(kmer_length),dtype=np.int64)

    
    with open(ref_file, 'r') as handle:
        for rec in SeqIO.parse(handle, 'fasta'):
            seq=''.join(rec.seq)
            #kmers, rev_kmers=seq2kmer(seq[:-3], kernel, tran_table, kmer_length)
            kmers, rev_kmers=seq2kmer(seq[:-13], kernel, tran_table, kmer_length)


    #print(len(kmers),len(np.unique(kmers)))
    unique_kmers, index, counts = np.unique(kmers,return_index=True,return_counts=True)

    ref_db_array_f=np.zeros([len(index),2], dtype=np.int64)
    ref_db_array_f[:,0]=unique_kmers
    ref_db_array_f[:,1]=index

    #print(len(rev_kmers),len(np.unique(rev_kmers)))
    unique_kmers, index, counts = np.unique(rev_kmers,return_index=True,return_counts=True)

    ref_db_array_r=np.zeros([len(index),2], dtype=np.int64)
    ref_db_array_r[:,0]=unique_kmers
    ref_db_array_r[:,1]=index
    
    return ref_db_array_f, ref_db_array_r, seq
    

if __name__ == '__main__':
    import os, psutil, time

    default_kmer_length=21
    example_ref_file='/home/liliandong/workspace/iSNV/DB/GCF_009858895.2_ASM985889v3_genomic.fna'

    T0=time.time()
    ref_db_array_f, ref_db_array_r, seq = build_ref_db(example_ref_file, default_kmer_length)
    T1=time.time()

    #np.save('temp/ref_db_array_f.npy', ref_db_array)
    #np.save('temp/ref_db_array_r.npy', ref_db_array)
    
    print('time using: ', T1-T0)
    print(u'RAM using %.4f GB' % (psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024) )