#!/bin/bash
#
#SBATCH --job-name=annotatedisorder_ba_be
#SBATCH --begin=now+8hours
#SBATCH --time=12:00:00

. /proj/bioinfo/users/x_oxasa/sp/sp/bin/activate

python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_bv  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_bv/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_bw  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_bw/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_bx  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_bx/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_by  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_by/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_bz  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_bz/