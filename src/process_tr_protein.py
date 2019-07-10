from __future__ import division
import argparse
import pandas as pd

"""
    ID                            MSA  begin  pvalue  l_effective  n  \
203  Q03S07  LHRVKRSILYSILD -KRPKRAKVYWFI-    497  0.0332           12  2   
204  Q03S07                      LTLG LTLM    431  0.0036            4  2   
205  Q03S07                      SILY SILD    503  0.0034            4  2   
206  Q03S07  LHRVKRSILYSILD -KRPKRAKVYWFI-    497  0.0332           12  2   

     n_effective      TRD  model  disordered_overlap   Entry  Length tr_type  
203            2  HHrepID  cpHMM                   0  Q03S07     675   small  
204            2  XSTREAM   None                   0  Q03S07     675   small  
205            2  XSTREAM   None                   0  Q03S07     675   small  
206            2  HHrepID   None                   0  Q03S07     675   small
"""
def transform_tr_to_residue(rows):
    residues = [0] * int(rows.iloc[0]['Length'])
    
    for index, row in rows.iterrows():
        tr_length=int(row['n_effective']) * int(row['l_effective'])
        for i in range(int(row['begin'])-1, int(row['begin'])+tr_length-1):
            residues[i] = 1
    return residues

def transform_row_to_residue(row):
    residues = [0] * int(row['Length'])
    tr_length=int(row['n_effective']) * int(row['l_effective'])
    for i in range(int(row['begin'])-1, int(row['begin'])+tr_length-1):
        residues[i] = 1
    return residues

def read_tr_annotations(input_file, swissprot_filename, file="tr_new.csv"):
    swissprot_df = pd.read_table(swissprot_filename)
    swissprot_df = swissprot_df[['Entry', 'Length']]
    tr_df = pd.read_csv(input_file, skipinitialspace=True)
    tr_df = pd.merge(tr_df, swissprot_df, right_on='Entry', left_on='ID', how='left')
    
    # add type column to tr_df

    tr_df['tr_type']=pd.cut(tr_df.l_effective, [0, 4, 15, 2000], right=False, labels=['micro', 'small', 'domain'])
    annotations_residues = []
    tr_types = []
    tr_ids = []
    tr_unitlengths = []
    # This is now per row and I need it to be per protein
    for protein_id in tr_df['ID'].unique():
        #annotations_residues.append(tr_df[tr_df.ID==protein_id].apply(transform_tr_to_residue, axis=1))
        protein_df = tr_df[tr_df.ID==protein_id]
        
        tr_type=protein_df.tr_type.unique()

        for t in tr_type:
            tr_types.append(t)
            tr_ids.append(protein_id)


            protein_t_df = protein_df[protein_df.tr_type==t]
            if len(protein_df)>1:
                annotations_residues.append(transform_tr_to_residue(protein_t_df))
            else:
               annotations_residues.append(transform_row_to_residue(protein_t_df))

            tr_unitlengths.append(protein_t_df.l_effective)
    
    result_df = pd.DataFrame({'tr': annotations_residues, 'type':tr_types, 'unit_length':tr_unitlengths, 'Entry':tr_ids})
    if file:
        result_df.to_csv('tr_new.csv')

    #return result_df




def read_commandline_arguments():
    parser = argparse.ArgumentParser(
        description='Process disorder')
    parser.add_argument('-t', '--tr_input', type=str,
                        help='The path to the tr file.')
    parser.add_argument('-s', '--swissprot', type=str,
                        help='Path to the swissprot file')
    
    return parser.parse_args()


if __name__ == '__main__':
    pars = read_commandline_arguments()
    
    tr_matrix = read_tr_annotations(pars.tr_input, pars.swissprot)
    