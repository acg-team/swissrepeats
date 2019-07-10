import requests
import urllib2
import json
import os
import glob

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

species = 'Saccharomyces+cerevisiae'
workingdir = "/proj/bioinfo/users/x_oxasa/sp/data/"


def get_disorder_info(disorder_data):
    disnum = 0
    for entry in disorder_data:
        if entry['ann'].lower() == 'd':
            disnum += (entry['end'] - entry['start']) + 1
    return disnum


def get_swissprot(swissprot_filename, outdir):
    swissprot_df = pd.read_table(swissprot_filename)
    uniprot_ids = swissprot_df['Entry'].values

    for protein_id in uniprot_ids:
        proteindir = outdir + protein_id + '/'

        if not os.path.exists(proteindir):
            try:
                url2 = ('http://mobidb.bio.unipd.it/ws/entries/' + protein_id
                        + '/disorder')
                response = requests.get(url2)
                print(url2)
                if response.status_code == requests.codes.ok:
                    disorder = response.json()
                    get_fasta(protein_id, disorder, outdir)
                else:
                    print(response.status_code)
            except requests.exceptions.RequestException as e:
                print(e)



def process_disorder(workingdir):
    outdir = workingdir + '/disorder_annotations'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    table_file = open(outdir + '/disorder.csv', 'w+')

    header_str = ('uniprotID,disnum,disnum_disprot,disnum_long,disnum_pdb,' +
                  'disnum_pdb_xray,disnum_pdb_nmr')

    all_predictors = ['iupl', 'iups', 'jronn', 'glo', 'espD', 'espN', 'espX',
                      'dis465', 'disHL', 'vsl']

    for pred in all_predictors:
        header_str += ',disnum_' + pred

    table_file.write(header_str + ',seqlength\n')

    for proteindir in glob.glob(workingdir + "/*"):
        protein_id = os.path.basename(proteindir)
        print(proteindir)
        disorder_file = proteindir + '/' + protein_id + '_disorder.json'
        if not os.path.exists(disorder_file):
            continue
        with open(disorder_file, 'r+') as f_dis:
            try:
                disorder = json.load(f_dis)
            except ValueError as e:
                print(e)
                print("no data for " + protein_id)
                continue

            # Look at full consensus
            disnum = 0
            if 'consensus' in disorder:
                full_consensus = disorder['consensus']['full']
                disnum = get_disorder_info(full_consensus)
                print(disorder['consensus'].keys())
                disprot = disorder['consensus']['disprot']
                disnum_disprot = get_disorder_info(disprot)

                exlong = disorder['consensus']['long']
                disnum_long = get_disorder_info(exlong)

                pdb = disorder['consensus']['pdb']
                disnum_pdb = get_disorder_info(pdb)

                pdb_xray = disorder['consensus']['pdb_xray']
                disnum_pdb_xray = get_disorder_info(pdb_xray)

                pdb_nmr = disorder['consensus']['pdb_nmr']
                disnum_pdb_nmr = get_disorder_info(pdb_nmr)

                disnum_predictors = {}
                for entry in disorder['predictors']:
                    disnum_predictors[entry['pred']] = get_disorder_info(entry['anns'])

            if (disnum > 0):
                # Read fasta file to know the sequnce length
                fasta_filename = proteindir + '/' + protein_id + '.fa'
                record = SeqIO.read(open(fasta_filename), 'fasta')
                seqlength = len(record.seq)

                protein_str = (protein_id +
                               ',' + str(disnum) +
                               ',' + str(disnum_disprot) +
                               ',' + str(disnum_long) +
                               ',' + str(disnum_pdb) +
                               ',' + str(disnum_pdb_xray) +
                               ',' + str(disnum_pdb_nmr))
                for pred in all_predictors:
                    if pred in disnum_predictors.keys():
                        protein_str += ',' + str(disnum_predictors[pred])
                    else:
                        protein_str += ',0'

                table_file.write(protein_str + ',' + str(seqlength) + '\n')
    table_file.close()


def process_disorder_coordinates(workingdir):
    outdir = workingdir + '/disorder_annotations'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    table_file = open(outdir + '/disorder_coordinates.csv', 'w+')
    header_str = 'uniprotID, start, end'
    table_file.write(header_str + '\n')

    for proteindir in glob.glob(workingdir + "/*"):
        protein_id = os.path.basename(proteindir)
        print(proteindir)
        disorder_file = proteindir + '/' + protein_id + '_disorder.json'
        if not os.path.exists(disorder_file):
            continue
        with open(disorder_file, 'r+') as f_dis:
            try:
                disorder = json.load(f_dis)
            except ValueError as e:
                print(e)
                print("no data for " + protein_id)
                continue

            # Look at full consensus
            if 'consensus' in disorder:
                full_consensus = disorder['consensus']['full']
                for entry in full_consensus:
                    if entry['ann'].lower() == 'd':
                        protein_str = (protein_id + "," + str(entry['start']) +
                                       ',' + str(entry['end']))
                        table_file.write(protein_str + '\n')
    table_file.close()



def process_disorder_regions(workingdir, swissprot_filename, method='consensus'):
    outdir = workingdir + '/disorder_annotations'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    swissprot_df = pd.read_table(swissprot_filename)
    uniprot_ids = swissprot_df['Entry'].values

    table_file = open(outdir + '/disorder_regions_'+ method + '.csv', 'w+')
    header_str = 'uniprotID, disorder, taxon'
    table_file.write(header_str + '\n')

    for proteindir in glob.glob(workingdir + "/*"):
        protein_id = os.path.basename(proteindir)
        print(proteindir)
        disorder_file = proteindir + '/' + protein_id + '_disorder.json'
        if not os.path.exists(disorder_file):
            continue
        with open(disorder_file, 'r+') as f_dis:
            try:
                disorder = json.load(f_dis)
            except ValueError as e:
                print(e)
                print("no data for " + protein_id)
                continue

            disorder_str = ''
            # Look at full consensus
            if 'consensus' in disorder:
                disorder_anns = []
                if method == 'consensus':
                   disorder_anns = disorder['consensus']['full']
                else:
                   predictors = disorder['predictors']
                   for entry in predictors:
                        print(entry['pred'])
                        if entry['pred'] == method:
                            disorder_anns = entry['anns']

                for entry in disorder_anns:
                        if entry['ann'].lower() == 'd':
                            disorder_str += (str(entry['start']) +
                                               ':' + str(entry['end']) + ';')
   
            swissprot_row = swissprot_df[swissprot_df['Entry'] == protein_id]
            protein_str = (protein_id + "," + disorder_str +
                               "," + swissprot_row['Taxonomic lineage (SUPERKINGDOM)'].iloc[0])
            print(protein_str)
            table_file.write(protein_str + '\n')
    table_file.close()


def get_fasta(proteinID, disorder, outdir):
    """
    Takes uniprot id, save a fasta file with a sequence
    and returns sequence length
    """

    try:
        url = "http://mobidb.bio.unipd.it/ws/entries/" + proteinID + "/uniprot"
        response = requests.get(url)
        print(url)
        if response.status_code == requests.codes.ok:
            data = response.json()
        else:
            print(response.status_code)
            return 0
    except requests.exceptions.RequestException as e:
        print(e)
        return 0

    protein = SeqRecord(Seq(data["sequence"]), id=proteinID, description="")
    proteindir = outdir + proteinID+"/"
    if not os.path.exists(proteindir):
        os.makedirs(proteindir)
    else:
        return len(data["sequence"])

    f_result = open(proteindir+proteinID+".fa", "w+")
    SeqIO.write(protein, f_result, "fasta")
    f_result.close()

    with open(proteindir+proteinID+"_disorder.json", "w+") as f_dis:
        json.dump(disorder, f_dis)

    return len(data["sequence"])


def annotationToString(entries, flag = 'd'):
    ids = []
    for entry in entries: 
        if entry['ann'].lower() == flag:
            for i in xrange(int(entry['start']), int(entry['end'])+1): 
                ids.append(i)
    return ids

def classifyDisorder():
     for proteindir in os.listdir(workingdir):
        with open(workingdir + proteindir+"/"+proteindir+"_disorder.json") as f_dis:
            disorder = json.load(f_dis) 
            
            disprot = annotationToString(disorder['consensus']['disprot'])
            pdb = annotationToString(disorder['consensus']['pdb']) #just disorder
            pdbs = annotationToString(disorder['consensus']['pdb'],'s') #just order
            
            pred = annotationToString(disorder['consensus']['predictors'])
            
            disseq = ""
     
            full_consensus = disorder['consensus']['full']
            for entry in full_consensus: 
                if entry['ann'].lower() == 'd':
                    for i in xrange(int(entry['start'])-1, int(entry['end'])): 
                        if ((i+1) in disprot): 
                            if ((i+1) in pdb):
                                disseq += '0'
                            elif ((i+1) in pdbs):
                                disseq += '1'
                            else:
                                disseq += '0'
                        else: 
                            if ((i+1) in pdb):
                                disseq += '4'
                            elif ((i+1) in pred):
                                disseq += '5'
                else:
                    for i in xrange(int(entry['start'])-1, int(entry['end'])): 
                        disseq += '1'
           
            fout = open(workingdir + proteindir+"/"+proteindir+"_disorder.fasta", "w+")       
            record = SeqRecord(id = proteindir,seq = Seq(disseq))

            SeqIO.write(record, fout, "fasta")
            fout.close()
