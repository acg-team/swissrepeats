#!/bin/bash
#
#SBATCH --job-name=annotatedisorder
#SBATCH --begin=now+3hours
#SBATCH --time=12:00:00

. /proj/bioinfo/users/x_oxasa/sp/sp/bin/activate

python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_ao  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_ao/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_ap  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_ap/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_aq  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_aq/