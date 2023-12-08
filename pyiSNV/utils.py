#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 22:03:37 2023

@author: lee
"""

import numpy as np
from hirola import HashTable

def encode(seq, kernel, tran_table):
    num_seq=seq.translate(tran_table)
    num_seq=np.array(list(num_seq),dtype=np.int64)

    #print(len(seq), len(kernel), len(tran_table))
    convolved=np.convolve(num_seq,kernel,mode='valid')
    
    assert len(convolved)==1
    return convolved[0]

def decode(kmer, kmer_length=21, return_base=True):
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

def build_tran_table_old():
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

def build_tran_table():
    tran_table={}
    for i in range(65,91):
        tran_table[chr(i)]='4'
    #for i in range(97,123):
    #    tran_table[chr(i)]='4'
    tran_table['A']='0'
    tran_table['C']='1'
    tran_table['G']='2'
    tran_table['T']='3'
    tran_table['\n']='4'
    return str.maketrans(tran_table)

def dict2table_old(the_dict, key_dtype=np.int64, value_dtype=np.int64, length=0):
    if length==0:
        length=int(len(the_dict)*2)
    table = HashTable(
        length,  # <--- Maximum size for the table 
        key_dtype,  # <--- NumPy dtype 
    )

    #value_array=np.zeros(len(the_dict),dtype=value_dtype)
    table.add(np.array(list(the_dict.keys())))
    value_array=np.array(list(the_dict.values()))
    #for key in the_dict.keys():
    #    idx=table.add(key_dtype(key))
    #    value_array[np.array(idx,dtype=value_dtype)]=the_dict[key]

    return table, value_array

def dict2table(the_dict, key_dtype=np.int64, value_dtype=np.int64, length=0, threshold=0):
    if length==0:
        length=int(len(the_dict)*2)
    table = HashTable(
        length,  # <--- Maximum size for the table 
        key_dtype,  # <--- NumPy dtype 
    )

    key_array=np.array(list(the_dict.keys()))
    value_array=np.array(list(the_dict.values()))
    if threshold > 0:
        key_array=key_array[value_array>=threshold]
        value_array=value_array[value_array>=threshold]
    #value_array=np.zeros(len(the_dict),dtype=value_dtype)
    table.add(key_array)
    #for key in the_dict.keys():
    #    idx=table.add(key_dtype(key))
    #    value_array[np.array(idx,dtype=value_dtype)]=the_dict[key]

    return table, value_array


tran_table=build_tran_table()