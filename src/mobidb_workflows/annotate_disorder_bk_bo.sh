#!/bin/bash
#
#SBATCH --job-name=annotatedisorder_ba_be
#SBATCH --begin=now+6hours
#SBATCH --time=12:00:00

. /proj/bioinfo/users/x_oxasa/sp/sp/bin/activate

python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_bk  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_bk/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_bl  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_bl/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_bm  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_bm/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_bn  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_bn/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_bo  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_bo/