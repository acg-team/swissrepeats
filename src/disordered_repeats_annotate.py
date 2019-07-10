from __future__ import division
import argparse
import csv
import os
from collections import defaultdict
import pandas as pd
import numpy as np
from Bio import SeqIO
from nucleotide_counts import read_regions
from repeat_disorder_overlap import find_overlap

def get_tr_characteristics(tr_row):
    id = tr_row['ID']
    begin = tr_row['begin']
    length =  tr_row['l_effective'] * tr_row['n']
    end = length + begin

    return id, begin, end, length

def compute_idr_percentage(tr_df, regions):
    idr_perc=[]
    for ind, row in tr_df.iterrows():
        id, begin, end, length = get_tr_characteristics(row)

        disordered_count = find_overlap(regions, id, begin, end)
        disorder_perc = disordered_count/length
       
        if disorder_perc > 1:
            print(id)
            print(disordered_count)
            print(length)
        
        idr_perc.append(disorder_perc)

    return idr_perc


def transform_fasta_to_dict(fasta):
    res_dict={}
    for record in fasta:
        #format 'sp|Q6GZX4|001R_FRG3G'
        uniprot_id=record.id.split('|')[1]
        seq = str(record.seq)
        res_dict[uniprot_id]=seq
    return res_dict

def compute_avg_topIDP(fasta_dict, tr_df):
    avg_topIDP=[]
    topIDP_scale = pd.read_table('data/scales_and_slopes.csv', index_col=0)
    print(topIDP_scale.loc['A','TopIDP'])

    for ind, row in tr_df.iterrows():
        id, begin, end, length = get_tr_characteristics(row)
        
        avg_topidp = np.NaN

        if id in fasta_dict:
            seq_region = fasta_dict[id][begin:end]
            sum_idp=0
            for aa in seq_region:
                if aa in topIDP_scale.index.values:
                    sum_idp+=topIDP_scale.loc[aa,'TopIDP']
            if len(seq_region)==0:
                avg_topidp = np.NaN
            else:
                avg_topidp = sum_idp/len(seq_region)
        
        avg_topIDP.append(avg_topidp)
    return avg_topIDP

def compute_aa_freqs(fasta_dict, tr_df):
    aminoacids = 'ABCDEFGHIKLMNOPQRSTUVWXYZ'
    aa_frequencies = {k:[] for k in aminoacids}
    print(aa_frequencies)
    for ind, row in tr_df.iterrows():
        id, begin, end, length = get_tr_characteristics(row)
        #aa_freqs=dict.fromkeys(aminoacids, 0)
        if id in fasta_dict:
            seq_region = fasta_dict[id][begin:end]
            for aa in aminoacids:
                aa_freq = seq_region.count(aa)/length
                aa_frequencies[aa].append(aa_freq)
    return aa_frequencies

def get_full_repeat_sequence(fasta_dict, tr_df):
    seqs = []
    for ind, row in tr_df.iterrows():
        id, begin, end, length = get_tr_characteristics(row)
        #aa_freqs=dict.fromkeys(aminoacids, 0)
        if id in fasta_dict:
            seq_region = fasta_dict[id][begin:end]
            seqs.append(seq_region)
    return seqs


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--repeats', required=True,
                        help='path to the repeat description file')
    parser.add_argument('-a', '--disorder', required=True,
                        help='path to the disoder file')
    parser.add_argument('-f', '--fasta', required=True,
                        help='path to the input fasta file')
   
    args = parser.parse_args()
    repeats_file = args.repeats
    disorder_file = args.disorder
    fasta_file = args.fasta
    regions, kingdoms = read_regions(disorder_file, 0)
    tr_df = pd.read_csv(repeats_file, index_col=False)
    tr_df = tr_df.dropna(axis=0, how='any')
    
    #idr_perc = compute_idr_percentage(tr_df, regions)
    #tr_df['idr_perc']=idr_perc
    #print(tr_df.head(10))

    fasta = SeqIO.parse(open(fasta_file, 'r'), 'fasta')
    fasta_dict=transform_fasta_to_dict(fasta)
    
    #avg_idp = compute_avg_topIDP(fasta_dict, tr_df)
    #tr_df['avg_idp']=avg_idp
    #print(tr_df.head(10))
    
    #aa_freqs = compute_aa_freqs(fasta_dict, tr_df)
    #aa_freqs_df = pd.DataFrame.from_dict(aa_freqs)
    #print(aa_freqs_df.head(10))
    #tr_df = pd.concat([tr_df, aa_freqs_df], axis=1)
    #print(tr_df.head(10))
    
    full_seqs = get_full_repeat_sequence(fasta_dict, tr_df)
    tr_df['full_seq']=full_seqs 

    for ind, row in tr_df.iterrows():
        print(">"+str(ind))
        print(row['full_seq'])
    #tr_df.to_csv('results/tr_annotations/tr_all_with_full_seq.csv', index=False)