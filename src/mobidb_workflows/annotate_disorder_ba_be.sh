#!/bin/bash
#
#SBATCH --job-name=annotatedisorder_ba_be
#SBATCH --begin=now+5hours
#SBATCH --time=12:00:00

. /proj/bioinfo/users/x_oxasa/sp/sp/bin/activate

python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_ba  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_ba/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_bb  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_bb/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_bc  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_bc/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_bd  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_bd/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_be  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_be/