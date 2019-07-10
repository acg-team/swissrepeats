####


## Create main overview Table 1

The main overview table 1 contains summary statistics on TR annotations, disordered TR annotations, and disorder annotations.
The script `src/create_tr_disorder_overlap_table.py` is used to compute these summary statistics.

It used `Sequence` and `Sequences` classes defined in `src/sequence_annotation_overlap.py`. This were designed to I/O multiple types of annotations,
and compute overlaps and summary statistics.
