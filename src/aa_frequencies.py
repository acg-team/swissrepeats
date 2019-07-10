import argparse
from Bio import SeqIO
import pandas as pd
import numpy as np

aminoacids = 'ABCDEFGHIKLMNOPQRSTUVWXYZ'


def get_fasta_data_handle(fasta_file):
    handle = SeqIO.parse(open(fasta_file, 'r'), 'fasta')
    return handle


def get_counts(fasta):
    counts = dict()
    for record in fasta:
        seq = str(record.seq)
        for aa in aminoacids:
            if aa in counts:
                counts[aa] += seq.count(aa)
            else:
                counts[aa] = seq.count(aa)
    return counts

def get_counts_region(fasta, regions):
    counts = dict()
    for record in fasta:
        #format 'sp|Q6GZX4|001R_FRG3G'
        uniprot_id=record.id.split('|')[1]
        seq = str(record.seq)
        if uniprot_id in regions:
            for region in regions[uniprot_id]:
                seq_region = seq[region[0]:region[1]]
                for aa in aminoacids:
                    if aa in counts:
                        counts[aa] += seq_region.count(aa)
                    else:
                        counts[aa] = seq_region.count(aa)
                print counts
    return counts

def get_repeats_regions(repeats_file):
    tr_df = pd.read_csv(repeats_file, skipinitialspace=True)
    regions=dict()
    for uniprot_id in tr_df.Entry.unique():
        print uniprot_id
        protein_rows = tr_df[tr_df.Entry==uniprot_id].tr
        for row in protein_rows:
            row = eval(row)
            start = np.argmax(row)
            row_rev = row[::-1]
            end = len(row_rev) - np.argmax(row_rev) - 1

            if uniprot_id in regions:
                regions[uniprot_id].append((start,end))
            else:
                regions[uniprot_id]=[(start,end)]

    return regions



def read(file):
    with open(file, 'rb') as fh:
        return pickle.load(fh)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', required=True,
                        help='path to the input fasta file')
    parser.add_argument('-a', '--annotations', required=True,
                        help='path to the annotations file')
    parser.add_argument('-o', '--output', required=True,
                        help='path to the output folder')
    args = parser.parse_args()
    fasta_file = args.fasta
    annotations_file = args.annotations

    fasta = get_fasta_data_handle(fasta_file)
    # All swissprot
    counts_dict=get_counts(fasta)
    counts_df=pd.DataFrame(counts_dict, index=[0])
    counts_df.to_csv("results/aa_freqs_all.csv")

    # All repeats
    #regions = get_repeats_regions('results/tr_annotations/tr_new.csv')
    #counts_dict=get_counts_region(fasta, regions)
    #counts_df=pd.DataFrame(counts_dict, index=[0])
    #counts_df.to_csv("results/aa_freqs_tr_domain.csv")

  

