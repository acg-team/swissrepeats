sleep 1h
python src/annotate_disorder.py -i /big/oxanas/swissrepeat/data/uniprot/split/swissprot_aw -o  /big/oxanas/swissrepeat/data/mobidb/split/swissprot_aw/ >> swissprot_aw_job.txt &
python src/annotate_disorder.py -i /big/oxanas/swissrepeat/data/uniprot/split/swissprot_ax -o  /big/oxanas/swissrepeat/data/mobidb/split/swissprot_ax/ >> swissprot_ax_job.txt &
python src/annotate_disorder.py -i /big/oxanas/swissrepeat/data/uniprot/split/swissprot_ay -o  /big/oxanas/swissrepeat/data/mobidb/split/swissprot_ay/ >> swissprot_ay_job.txt &
python src/annotate_disorder.py -i /big/oxanas/swissrepeat/data/uniprot/split/swissprot_az -o  /big/oxanas/swissrepeat/data/mobidb/split/swissprot_az/ >> swissprot_az_job.txt
wait
python src/annotate_disorder.py -i /big/oxanas/swissrepeat/data/uniprot/split/swissprot_ca -o  /big/oxanas/swissrepeat/data/mobidb/split/swissprot_ca/ >> swissprot_ca_job.txt &
python src/annotate_disorder.py -i /big/oxanas/swissrepeat/data/uniprot/split/swissprot_cb -o  /big/oxanas/swissrepeat/data/mobidb/split/swissprot_cb/ >> swissprot_cb_job.txt &
python src/annotate_disorder.py -i /big/oxanas/swissrepeat/data/uniprot/split/swissprot_cc -o  /big/oxanas/swissrepeat/data/mobidb/split/swissprot_cc/ >> swissprot_cc_job.txt &
python src/annotate_disorder.py -i /big/oxanas/swissrepeat/data/uniprot/split/swissprot_cd -o  /big/oxanas/swissrepeat/data/mobidb/split/swissprot_cd/ >> swissprot_cd_job.txt
wait