## Swissprot

#### Retrieve the sequence data.
    cd /home/eschaper/archive/tral/sp_0/data/seq
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2015_11/knowledgebase/uniprot_sprot-only2015_11.tar.gz
    tar -xzf uniprot_sprot-only2015_11.tar.gz
    mv uniprot_sprot-only2015_11/uniprot_sprot.fasta.gz .
    gunzip uniprot_sprot.fasta.gz

### optional:
    rm -rf uniprot_sprot-only2015_11*

Split the sequence files:

	wget ftp://saf.bio.caltech.edu/pub/software/molbio/fastasplitn.c
	gcc -Wall -std=c99 -pedantic -o fastasplitn fastasplitn.c
	cd /scratch/cluster/monthly/eschaper/tral/sp/data/seq/split
	fastasplitn -in /scratch/cluster/monthly/eschaper/tral/sp/data/seq/uniprot_sprot.fasta -n 100


#### Retrieve PFAM models
Download Pfam-A

	wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz


Turn each HMM into a .pickle

	 python3 /home/eschaper/Python_projects/tral/tral/examples/workflow/tandem_repeat_annotation_scripts.py file_preparation --hmm_annotation_raw /scratch/cluster/monthly/eschaper/tral/sp/data/PFAM_annotations/uniprot_sprot_pfam.pickle --hmm_annotation /scratch/cluster/monthly/eschaper/tral/sp/data/PFAM_annotations/split --hmm_raw /scratch/cluster/monthly/eschaper/tral/sp/data/PFAM/Pfam-A.hmm --hmm /scratch/cluster/monthly/eschaper/tral/sp/data/PFAM/split


This could produce 16230 files (however, in the current version, sequence profile models with l>200 were ignored.):

	[eschaper@devfrt01 PFAM]$ cat Pfam-A.hmm | grep "HMMER3" | wc -l
	16230


#### Retrieve PFAM annotations
Go to http://www.uniprot.org/uniprot/?query=reviewed%3Ayes&sort=score
and download data in tab separated format. Then, copy the data to the cluster.

	scp uniprot-reviewed%3Ayes.tab.gz eschaper@prd.vital-it.ch:/scratch/cluster/monthly/eschaper/tral/sp/data/PFAM uniprot-reviewed%3Ayes.tab.gz


Remove those PFAM annotations relating to large PFAM models (l>200):
	
	python3 filter_pfam_annotations.py /scratch/cluster/monthly/eschaper/tral/sp/data/PFAM/uniprot-reviewed%3Ayes.tab /scratch/cluster/monthly/eschaper/tral/sp/data/PFAM/uniprot-reviewed%3Ayes_filtered.tab /scratch/cluster/monthly/eschaper/tral/sp/data/PFAM/split /scratch/cluster/monthly/eschaper/tral/sp/data/PFAM/Pfam_models_20151209_smaller_200chars.txt


Create `sp_pfam_annotations_20151208_smaller_200chars.pickle` from `Pfam_models_20151209_smaller_200chars.txt` using the method `read_pfam_uniprot` in `tandem_repeat_annotation_scripts.py`.


#### Test workflow manually

	export MYTRAL=/home/eschaper/Python_projects/tral/tral
	export MYGROUND=/scratch/cluster/monthly/eschaper/tral/sp
	python3 $MYTRAL/examples/workflow/tandem_repeat_annotation_scripts.py workflow -i $MYGROUND/data/seq/split/frag100 -o $MYGROUND/results/split/frag100.pickle -os $MYGROUND/results/serialized/frag100.tsv -f tsv -t 600  --hmm_annotation $MYGROUND/data/PFAM/sp_pfam_annotations_20151208.pickle --hmm $MYGROUND/data/PFAM/split

### GC3PIE: Local installation
Login to LSF cluster. Then:

	wget https://raw.githubusercontent.com/uzh/gc3pie/master/scripts/install.py
	
	. /home/eschaper/software/lib/gc3pie/bin/activate

	
#### Set up the local GC3PIE workflow config

Save to `tandem_repeat_annotation_workflow_sp.ini `:


	required_memory = 5000
	[file_preparation]
	    activated = False
	    script = python3 /home/eschaper/Python_projects/tral/tral/examples/workflow/tandem_repeat_annotation_scripts.py file_preparation
	    input = /scratch/cluster/monthly/eschaper/tral/sp/data/uniprot_PRDM_annotation.tsv
	    output = /scratch/cluster/monthly/eschaper/tral/sp/data/uniprot_PRDM_annotation.pickle
	    extra = --hmm_annotation_raw /scratch/cluster/monthly/eschaper/tral/sp/data/PFAM_annotations/uniprot_sprot_pfam.pickle --hmm_annotation /scratch/cluster/monthly/eschaper/tral/sp/data/PFAM_annotations/split --hmm_raw /scratch/cluster/monthly/eschaper/tral/trembl/data/PFAM/PFAM-A.hmm --hmm /scratch/cluster/monthly/eschaper/tral/trembl/data/PFAM/split
	    logdir = /scratch/cluster/monthly/eschaper/tral/sp/log/split_sequence_file
	    stdout = stdout.log
	    stderr = stderr.log
	[sequencewise_parallel_flow]
	    input = /scratch/cluster/monthly/eschaper/tral/sp/data/seq/split
	[annotate_tandem_repeats]
	    activated = True
	    script = python3 /home/eschaper/Python_projects/tral/tral/examples/workflow/tandem_repeat_annotation_scripts.py workflow
	    input = /scratch/cluster/monthly/eschaper/tral/sp/data/seq/split/$N
	    output = /scratch/cluster/monthly/eschaper/tral/sp/results/split/$N.pickle
	    extra = -os /scratch/cluster/monthly/eschaper/tral/sp/results/serialized/$N.tsv -f tsv -t 86400 --hmm_annotation /scratch/cluster/monthly/eschaper/tral/sp/data/PFAM/sp_pfam_annotations_20151208_smaller_200chars.pickle --hmm /scratch/cluster/monthly/eschaper/tral/sp/data/PFAM/split
	    logdir = /scratch/cluster/monthly/eschaper/tral/sp/log/annotate_tandem_repeats
	    stdout = stdout.log
	    stderr = stderr.log
	    required_memory = 80000	
	
	
#### Run the GC3PIE workflow

	source ~/software/lib/gc3pie/bin/activate
	export MYCODE=/home/eschaper/workflow
	export HOST=vital-it
	// You can add -u sqlite:////path/to/<session_name>.db
	export SESSION=dr100
	python2 $MYTRAL/examples/workflow/tandem_repeat_annotation_workflow.py -w 60minutes -r $HOST -J 500  -s $SESSION -C 2 -vvvv -conf $MYCODE/tandem_repeat_annotation_workflow_sp.ini -u sqlite:///$MYCODE/$SESSION.db


#### Download additional annotations from Swiss-Prot

http://www.uniprot.org/uniprot/?query=reviewed%3Ayes&sort=score&columns=id%2Cprotein%20names%2Clength%2Clineage(SUPERKINGDOM)%2Clineage(KINGDOM)%2Clineage(ORDER)%2Clineage(CLASS)%2Clineage(FAMILY)%2Clineage(SPECIES)%2Clineage-id(SPECIES)%2Cdatabase(DisProt)%2Cdatabase(MobiDB)%2Cvirus%20hosts%2Cdatabase(OrthoDB)%2Cdatabase(OMA)


###  To tandem repeat annotations, add how many amino acids in the tandem repeat region are annotated as disordered


	python /proj/bioinfo/users/x_oxasa/sp/src/repeat_disorder_overlap.py -r /proj/bioinfo/users/x_oxasa/sp/results/tr_annotations/concatenated.csv -a /proj/bioinfo/users/x_oxasa/sp/results/disorder_annotations/mobidb_regions.csv -o /proj/bioinfo/users/x_oxasa/sp/results/tr_annotations/trs_annotated_with_disorder


