#!/bin/bash
#
#SBATCH --job-name=annotatedisorder_ba_be
#SBATCH --begin=now+7hours
#SBATCH --time=12:00:00

. /proj/bioinfo/users/x_oxasa/sp/sp/bin/activate

python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_bp  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_bp/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_br  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_br/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_bs  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_bs/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_bt  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_bt/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_bu  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_bu/