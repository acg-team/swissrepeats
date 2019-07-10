#!/bin/bash
#
#SBATCH --job-name=annotatedisorder_j_n
#SBATCH --begin=now+2hours
#SBATCH --time=12:00:00

. /proj/bioinfo/users/x_oxasa/sp/sp/bin/activate

python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_aj  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_aj/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_ak  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_ak/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_al  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_al/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_am  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_am/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_an  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_an/