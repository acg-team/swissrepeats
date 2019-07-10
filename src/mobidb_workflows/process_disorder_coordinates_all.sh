#!/bin/bash
#
#SBATCH --job-name=process_disorder
#SBATCH --time=12:00:00
#SBATCH --ntasks=30

. /proj/bioinfo/users/x_oxasa/sp/sp/bin/activate

MOBIDB=/proj/bioinfo/users/x_oxasa/sp/data/mobidb/split/*

for file in $MOBIDB 
do
   (python /proj/bioinfo/users/x_oxasa/sp/src/process_disorder.py -t c -i "$file") &
  
done
wait

