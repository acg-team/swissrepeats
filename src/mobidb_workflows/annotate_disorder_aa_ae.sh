#!/bin/bash
#
#SBATCH --job-name=annotatedisorder
#SBATCH --begin=now+1hour
#SBATCH --time=12:00:00

. /proj/bioinfo/users/x_oxasa/sp/sp/bin/activate

python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_aa  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_aa/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_ab  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_ab/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_ac  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_ac/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_ad  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_ad/
python /proj/bioinfo/users/x_oxasa/sp/src/annotate_disorder.py -i /proj/bioinfo/users/x_oxasa/sp/data/uniprot/split/swissprot_ae  -o /proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/swissprot_ae/