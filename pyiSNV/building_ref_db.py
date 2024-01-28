# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:58:54 2022

@author: lilab
"""

import numpy as np
from Bio import SeqIO

from pyiSNV.utils import tran_table

class GenomeMapBuilder:
    def __init__(self, config):
        self.kmer_length = config.kmer_length
    
    def build_ref_db(self, ref_file):
        
        kmer_length = self.kmer_length
        
        if kmer_length%2 == 0:
            print('ERROR: setting kmer length as a even number is currently not supported')
            kmer_length = kmer_length + 1
            print('Revised kmer length:', kmer_length)  
            
        if kmer_length > 29:
            print('ERROR: kmer length larger than 29 is currently not supported')
            kmer_length = 29
            print('Revised kmer length:', kmer_length)
        elif kmer_length < 1:
            print('ERROR: invlid kmer length')
            kmer_length = 21
            print('Revised kmer length:', kmer_length)

        kernel=4**np.array(range(kmer_length),dtype=np.int64)

        
        with open(ref_file, 'r') as handle:
            for rec in SeqIO.parse(handle, 'fasta'):
                seq=''.join(rec.seq)
                
                seq = remove_repeat(seq)
                #kmers, rev_kmers=seq2kmer(seq[:-3], kernel, tran_table, kmer_length)
                kmers, rev_kmers=seq2kmer(seq, kernel, tran_table, kmer_length)
                if len(kmers) != len(seq) - kmer_length + 1:
                    print('ERROR: Depulicate kmer exist, please increase kmer length')
                assert len(kmers) == len(np.unique(kmers))


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
        
        ref_db_array=np.concatenate([ref_db_array_f, ref_db_array_r])

        self._ref_kmer_dict = dict(ref_db_array)
        self._ref_kmer_f_dict = dict(ref_db_array_f)
        self._ref_kmer_r_dict = dict(ref_db_array_r)
        
        self._genome = seq
        
        return True
    
    def get_ref_kmer(self, keyword='all'):

        if keyword == 'forward':
            return self._ref_kmer_f_dict
        elif keyword == 'reverse':
            return self._ref_kmer_r_dict
        else:
            return self._ref_kmer_dict
    
    def get_genome_seq(self):
        return self._genome
    
def remove_repeat(seq):
    for i in range(len(seq)):
        if seq[-i - 1] == 'A' or seq[-i - 1] == 'a':
            continue
        else:
            #print(i)
            return seq[:-i]

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


def build_ref_db_bk(ref_file, kmer_length=21):

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
    
