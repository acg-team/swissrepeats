from __future__ import division
import argparse
import pandas as pd
import numpy as np
import os
import math

"""
Step1: 

Calculate matrix
    method1 method2 ...
res1  1      1
res2  1      1
res3  1      0
...
"""

def transform_to_residue(row):
    residues = [0] * int(row['Length'])
    if not (pd.isnull(row['disorder'])):
        for ann in row['disorder'].split(';'):
            ann = ann.split(':')
            if len(ann) > 1:
                for i in range(int(ann[0])-1,int(ann[1])):
                    if i>=row['Length']:
                        residues=[]
                        continue
                    else:
                        residues[i] = 1
    return residues


def transform_tr_to_residue(row):
    residues = [0] * int(row['Length'])
    tr_length=int(row['n_effective']) * int(row['l_effective'])
    for i in range(int(row['begin'])-1, int(row['begin'])+tr_length-1):
        residues[i] = 1
    return residues

def read_disorder_annotations(input_folder, swissprot_filename, method):
    """
      Transform: 
      P25670,1:6;8:15;46:56;59:68;,Eukaryota
      P63115,1:4;6:6;13:22;424:430;489:494;508:530;,Eukaryota
      Into:
      1
      2
      3
      4
      5
      6
    """
    methods=[method]
    #methods=['iupl','iups','jronn', 'dis465', 'disHL', 'espX', 'vsl', 'espD', 'espN']

    swissprot_df = pd.read_table(swissprot_filename)
    swissprot_df = swissprot_df[['Entry', 'Length']]

    print swissprot_df.head(5)

    all_annotations = {}
    for method in methods:
        method_df = pd.read_csv(input_folder + 'mobidb_regions_' + method + '.csv', skipinitialspace=True)
        print len(method_df)
        all_annotations[method] = method_df['disorder'].values

    all_annotations['Entry'] = method_df['uniprotID'].values

    annotations_df = pd.DataFrame.from_dict(all_annotations) 
    
    # Do a left merge with swissprot based on Entry column 
    annotations_df = pd.merge(annotations_df, swissprot_df, on='Entry', how='left')

    print annotations_df.columns.values

    annotations_residues = {}
    for method in methods:
       if os.path.isfile(method+'.csv'):
            annotations_residues[method] = pd.read_csv(method+'.csv', names=["index", "ann"]).ann
            print annotations_residues[method].head(5)
            continue

       temp_df = pd.DataFrame({'disorder': annotations_df[method], 'Length':annotations_df['Length']})
       if method in annotations_residues:
          annotations_residues[method].append(temp_df.apply(transform_to_residue, axis=1))
       else:
          annotations_residues[method]=temp_df.apply(transform_to_residue, axis=1)
       print annotations_residues[method].head(5)
       annotations_residues[method].to_csv(method+'.csv')
 
    result_df = pd.DataFrame.from_dict(annotations_residues)
    result_df['Entry'] = all_annotations['Entry']
    print result_df.head(5)
    return result_df

def read_tr_annotations(input_file, swissprot_filename):
    swissprot_df = pd.read_table(swissprot_filename)
    swissprot_df = swissprot_df[['Entry', 'Length']]
    tr_df = pd.read_csv(input_file, skipinitialspace=True)
    tr_df = pd.merge(tr_df, swissprot_df, right_on='Entry', left_on='ID', how='left')
    
    # add type column to tr_df

    tr_df['tr_type']=pd.cut(tr_df.l_effective, [0, 4, 15, 2000], right=False, labels=['micro', 'small', 'domain'])
    annotations_residues = {}

    for t in set(tr_df.tr_type):
        annotations_residues['tr_'+t]=tr_df[tr_df.tr_type==t].apply(transform_tr_to_residue, axis=1)
   
    result_df = pd.DataFrame.from_dict(annotations_residues)
    result_df['Entry'] = tr_df['ID'] 

    return result_df


def mathews_row(row):
    from sklearn.metrics import matthews_corrcoef

    if type(row['method1']) is str:
        row['method1']=eval(row['method1'])
    if type(row['method2']) is str:
        row['method2']=eval(row['method2'])

    if (np.isnan(row['method1']).all()) or (np.isnan(row['method2']).all()):
       mcc_residues = 0
    elif ((len(row['method1'])<1) or (len(row['method2'])<1)):
        mcc_residues = 0
    else:
        mcc_residues = matthews_corrcoef(row['method1'], row['method2'])  
    mcc_protein = np.mean(mcc_residues)
    return mcc_protein

def calculate_mcc(ann_df):
    methods = ann_df.columns.values
    all_disorder = ['iupl','iups', 'jronn', 'dis465', 'disHL', 'espX', 'vsl', 'espD', 'espN']
    all_tr = ['tr_domain', 'tr_small', 'tr_micro']
    itemindex = np.where(methods=='Entry')
    methods=np.delete(methods, itemindex)

    mcc_matrix = np.zeros(shape=(len(methods),len(methods)), dtype=np.float32)
    mcc_df = pd.DataFrame(mcc_matrix, columns=methods)
    mcc_df.set_index(mcc_df.columns.values, inplace=True)
    
    for i in methods:
        for j in methods:
            if ((i==j) or ((i in all_tr) and (j in all_tr)) or ((i in all_disorder) and (j in all_disorder))): # this needed when we only want to correlate to tr
                mcc_df[i][j] = 1
                continue
            if mcc_df[i][j]!=0 or mcc_df[j][i]!=0:
                continue
            temp_df = pd.DataFrame({'method1': ann_df[i], 'method2':ann_df[j]})
            mcc_proteins=temp_df.apply(mathews_row, axis=1)
            mcc_total=mcc_proteins.sum()
            print mcc_total

            mcc_mean = mcc_total/len(ann_df)
            print str(len(ann_df))
            print mcc_mean
            mcc_df[i][j] = mcc_mean
    print mcc_df

def compute_mcc_numbers(row):
    #   calculate 4 numbers
    #   method1_and_method2=number of residues where method1==method2==1
    #   method1_not_method2=number of residues where method1==1 and method2==0
    #   not_method1_not_method2=number of residues where method1==method2==0
    #   not_method1_method2=number of residues where method1==0 and method2==1

    method1_and_method2=0
    method1_not_method2=0
    not_method1_not_method2=0
    not_method1_method2=0

    if type(row['method1']) is str:
        row['method1']=eval(row['method1'])
    if type(row['method2']) is str:
        row['method2']=eval(row['method2'])

    if (np.isnan(row['method1']).all()):
        row['method1'] = [0] * row['Length']  
    if (np.isnan(row['method2']).all()):
        row['method2'] = [0] * row['Length']

    # Speeding up tips: 
    #   calculate the sum, if it equals zero for both, then not_method1_not_method2 increment
    # compare two binary arrays: put it in a matrix? calculate contingency matrix for two arrays

    """
    for ind in range(0, len(row['method1'])):
        i = row['method1'][ind]
        j = row['method2'][ind]
        if i==1: 
            if j==1:   # 1 1
                method1_and_method2+=1
            else:  # 1 0 
                method1_not_method2+=1
        elif j==1:  # 0 1
            not_method1_method2+=1
        else:  # 0 0 
            not_method1_not_method2+=1
    """
    from sklearn.metrics import confusion_matrix
    #print row['method1']
    #print row['method2']

    confusion = confusion_matrix(row['method1'], row['method2'], [0, 1])

    #return [method1_and_method2, method1_not_method2, not_method1_not_method2, not_method1_method2]
    return np.asarray(confusion).reshape(-1).tolist()

def calculate_mcc_pairwise_matrix(ann_df, method1, method2):
    mcc_matrix = np.zeros(shape=(2,2), dtype=np.float32)
    mcc_df = pd.DataFrame(mcc_matrix, columns=[method1, method2])
    mcc_df.set_index(mcc_df.columns.values, inplace=True)
    
    temp_df = pd.DataFrame({'method1': ann_df[method1], 'method2':ann_df[method2], 'Length':ann_df['Length'] })
    
    print temp_df.head(5)
    print len(temp_df)
    # For the entire temp_df
   
    mcc_proteins=temp_df.apply(compute_mcc_numbers, axis=1)

    mcc_proteins = mcc_proteins.tolist()
    print "RESULTS"
    column_sums = [sum([row[i] for row in mcc_proteins]) for i in range(0,len(mcc_proteins[0]))]

    print "nothing, tr no disorder,  disorder no tr, tr+disorder"
    print column_sums
    print "% TR in IDR"
    print column_sums[3]/(column_sums[3]+column_sums[1])*100
    total_residues = sum(column_sums)
    print "% IDR"
    print (column_sums[2]+column_sums[3])/total_residues*100
    print "MCC"
    print ((column_sums[3]*column_sums[0])-(column_sums[1]*column_sums[2]))/math.sqrt((column_sums[3]+column_sums[1])*(column_sums[3]+column_sums[2])*(column_sums[0]+column_sums[1])*(column_sums[0]+column_sums[2]))

def read_commandline_arguments():
    parser = argparse.ArgumentParser(
        description='Process disorder')
    parser.add_argument('-i', '--input', type=str,
                        help='The path to the mobidb folder.')
    parser.add_argument('-t', '--tr_input', type=str,
                        help='The path to the tr file.')
    parser.add_argument('-s', '--swissprot', type=str,
                        help='Path to the swissprot file')
    parser.add_argument('--methodT', type=str,
                        help='TR method to correlate to')
    parser.add_argument('--methodD', type=str,
                        help='Disorder method to correlate to')
    parser.add_argument('--filter', type=str,
                        help='Optional filtering taxon: Eukaryota, Bacteria, Viruses')
    
    pars = vars(parser.parse_args())
    pars = {key: value for key, value in pars.iteritems() if value is not None}
    return pars

def read_tr_annotations_ready(input_file, method):
    tr_df = pd.read_csv(input_file, index_col=False)
    print "Total tr annotations"
    print len(tr_df)

    tr_filtered = tr_df[tr_df.type==method][['tr', 'Entry']]
    tr_filtered.columns=['tr_'+method, 'Entry']
    return tr_filtered

if __name__ == '__main__':
    pars = read_commandline_arguments()
    method1=pars['methodT']
    method2=pars['methodD']
    print method1
    print method2

    # filter only disorder matrix by taxon 
    disorder_matrix = read_disorder_annotations(pars['input'], pars['swissprot'], method2)
    disorder_matrix = disorder_matrix[[method2, 'Entry']].drop_duplicates()
   
    #tr_matrix = read_tr_annotations(pars['tr_input'], pars['swissprot'])
    tr_matrix = read_tr_annotations_ready(pars['tr_input'], method1)
    print "Length of tr matrix"
    print len(tr_matrix)
    print "Unique tr protein annotations"
    print len(tr_matrix.Entry.unique())
    print "Length of disorder matrix"
    print len(disorder_matrix.dropna(how='any'))
    print "Unique disorder protein"
    print len(disorder_matrix.drop_duplicates())
    
    swissprot_df = pd.read_table(pars['swissprot'])

    if 'filter' in pars:
        swissprot_df = swissprot_df[swissprot_df['Taxonomic lineage (SUPERKINGDOM)']==pars['filter']]

    swissprot_df = swissprot_df[['Entry', 'Length']]
    print "All swissprot"
    print len(swissprot_df)
    annotations_matrix = pd.merge(swissprot_df, disorder_matrix, how='left', left_index=True, on='Entry')
    print len(annotations_matrix)
    annotations_matrix = pd.merge(annotations_matrix, tr_matrix, how='left', left_index=True, on='Entry')  # need to merge these two matrixes to swissprot to preserve all seqs!
    print len(annotations_matrix)


    # pairwise methods matrix, all residues
    calculate_mcc_pairwise_matrix(annotations_matrix, method2, 'tr_'+method1)
    