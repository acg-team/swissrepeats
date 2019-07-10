import pandas as pd
import numpy as np
import argparse, itertools

def generate_fake_alignment():
	aa_freqs = pd.read_table('results/subjmeans.csv', index_col=0, sep=',')

	N = 1000
	alignment=[]
	aminoacids=aa_freqs.columns

	for row in aa_freqs.iterrows():
		index, data = row
		probs = data.tolist()
		probs /= sum(probs)

		col=np.random.choice(aminoacids, N, p=probs)

		alignment.append(col)

	msa = ([list(x) for x in zip(*alignment)])
	for row in msa:
		print(''.join(row))

def clstr_filt(cluster_dic, minimum):

  # filter the cluster in the cluster_dic based
  # on the minimum cluster size

  # grab a list of the cluster_dic keys
  # and create an empty list for the filtered sequence ids
  clusters, filtered_dic = list(cluster_dic), {}

  # parse through the keys and check the cluster size,
  # if the size is > minimum size, add the cluster to the
  # filtered dic
  for cluster in clusters:
    if len(cluster_dic[cluster]) >= minimum:
      for sequence in cluster_dic[cluster]:
        filtered_dic[sequence] = cluster

  # return the filtered cluster_dic
  return filtered_dic

def read_clstr(cluster_filename):

  # parse through the .clstr file and create a dictionary
  # with the sequences per cluster

  # open the cluster file and set the output dictionary
  cluster_file, cluster_dic = open(cluster_filename), {}

  # parse through the cluster file and store the cluster name + sequences in the dictionary
  cluster_groups = (x[1] for x in itertools.groupby(cluster_file, key=lambda line: line[0] == '>'))
  for cluster in cluster_groups:
    name = cluster.__next__().strip()
    seqs = [seq.split('>')[1].split('...')[0] for seq in cluster_groups.__next__()]
    #print(name)
    #print(seqs)
    cluster_dic[name] = seqs

  # return the cluster dictionary
  return cluster_dic


def get_seqs_from_clusters(cluster_filename, minimum):

	# obtain a dictionary with the clusters and sequences
	cluster_dic = read_clstr(cluster_filename)
	filtered_dic = clstr_filt(cluster_dic, minimum)

	x = list(itertools.islice(filtered_dic.items(), 0, 4))
	print(x)



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Filter the output from CD-hit based on the minimum number of read per cluster.\nThe filtered output fasta file produced has the same name as the input file with  _min_[minimum size from -c argument].fasta attachted to the name.')

	parser.add_argument('-f', '--fasta', metavar='.fasta file', dest='fasta', type=str,
	      help='The .fasta file containing the clusters produced by CD-hit.')
	parser.add_argument('-c', '--cluster', metavar='.clstr file', dest='cluster', type=str,
	      help='The .clstr file producec by CD-hit that contains the cluster information.')
	parser.add_argument('-m', '--minimum', metavar='minimum size', dest='minimum', type=int,
	      help='The minimum cluster size.')

	args = parser.parse_args()
	get_seqs_from_clusters(args.cluster, args.minimum)