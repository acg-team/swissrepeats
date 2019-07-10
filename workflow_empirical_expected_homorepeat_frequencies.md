##  Empirical and estimated homorepeat frequencies.

### Calculate the frequencies of amino acids in annotated regions for a sequence base

Calculate the frequencies of amino acids in disordered regions in Swiss-Prot:

    export MYDIR=/proj/bioinfo/users/x_oxasa/sp
    # Repeat for every type of disorderd regions
    export MYTAG=mobidb_regions_dis465
    export MYTAG=mobidb_regions_disHL
    export MYTAG=mobidb_regions_espD
    export MYTAG=mobidb_regions_espN
    export MYTAG=mobidb_regions_espX
    export MYTAG=mobidb_regions_inverted
    export MYTAG=mobidb_regions_iupl
    export MYTAG=mobidb_regions_iups
    export MYTAG=mobidb_regions_jronn
    export MYTAG=mobidb_regions_vsl
    export MYTAG=mobidb_regions

    python $MYDIR/src/nucleotide_counts.py -f $MYDIR/data/seq/uniprot_sprot.fasta -a $MYDIR/results/disorder_annotations/$MYTAG.csv -o $MYDIR/results/swiss_prot_aa_frequencies/$MYTAG


	
Use `create_regions_file_from_annotations_file` in `nucleotide_counts.py` to create an input file for `nucleotide_counts.py` allowing to calculate the frequency of amino acids in entire fasta files.

	create_regions_file_from_annotations_file("$MYDIR/data/uniprot/swissprot_annotations.tab", "$MYDIR/data/uniprot/swissprot_annotations_for_aa_counts.csv")
	
	
With this, calculate the frequencies of amino acids in the entire Swiss-Prot:

	python $MYDIR/src/nucleotide_counts.py -f $MYDIR/data/seq/uniprot_sprot.fasta -a $MYDIR/data/uniprot/swissprot_annotations_for_aa_counts.csv -o $MYDIR/results/swiss_prot_aa_frequencies/swissprot_kingdomwise

#### Calculate amino acids frequencies for ordered regions.
Also calculate the frequencies of amino acids in all regions which are not annotated as disordered. We assume these to be inverted compared to ordered regions.
To create the region annotation file for ordered regions; in the `src` directory, start Python3 and run:

    import os
    import nucleotide_counts
    path_trunk = "/path/to/swissrepeat"
    tags = ["mobidb_regions_dis465", "mobidb_regions_disHL", "mobidb_regions_espD", "mobidb_regions_espN", "mobidb_regions_espX", "mobidb_regions_iupl", "mobidb_regions_iups", "mobidb_regions_jronn", "mobidb_regions_vsl", "mobidb_regions"]
    tag = tags[4]
    # Load region annotations
    for tag in tags[4:]:
        regions_file = os.path.join(path_trunk, "results", "disorder_annotations", "{}.csv".format(tag))
        inverted_regions, kingdoms = nucleotide_counts.read_regions(regions_file, min_region_length=1)
        annotation_input_file = os.path.join(path_trunk, "data/swissprot_annotations.tsv")
        regions_result_file = os.path.join(path_trunk, "results", "disorder_annotations", "{}_inverted.csv".format(tag))
        delimiter_input="\t"
        delimiter_result=","
        nucleotide_counts.create_regions_file_from_annotations_file(annotation_input_file,regions_result_file,delimiter_input,delimiter_result,inverted_regions)

Next, calculate the frequencies of amino acids in disordered regions in Swiss-Prot:

    export MYTAG=mobidb_regions_dis465_inverted
    export MYTAG=mobidb_regions_disHL_inverted
    export MYTAG=mobidb_regions_espD_inverted
    export MYTAG=mobidb_regions_espN_inverted
    export MYTAG=mobidb_regions_espX_inverted
    export MYTAG=mobidb_regions_iupl_inverted
    export MYTAG=mobidb_regions_iups_inverted
    export MYTAG=mobidb_regions_jronn_inverted
    export MYTAG=mobidb_regions_vsl_inverted
    export MYTAG=mobidb_regions_inverted

    python $MYDIR/src/nucleotide_counts.py -f $MYDIR/data/seq/uniprot_sprot.fasta -a $MYDIR/results/disorder_annotations/$MYTAG.csv -o $MYDIR/results/swiss_prot_aa_frequencies/$MYTAG


	

### Calculate empirical homorepeat frequencies.


### Calculate expected homorepeat frequencies and connect empirical and expected data.

Requirements: 
- Annotations of empirical homorepeat counts
- Amino acid frequencies 


	# Ordered regions
	python3 homorepeat_counts.py -r results/disorder_annotations/mobidb_regions_inverted.csv -t swissprot_disorder_inverse_kingdomwise -p path/to/swissrepeat
	
	# Disordered regions
	python3 homorepeat_counts.py -r results/disorder_annotations/mobidb_regions.csv -t swissprot_disorder_kingdomwise -p path/to/swissrepeat
