#!/usr/bin/env bash
cd app
for FILE in 1 1_to_3 3_to_5
do
 echo ${FILE}
 bowtie2 -a --local --mp 2,2 --rdg 10,2 --rfg 10,2 -f -x ./sample_data/reference/reference -U ./sample_data/reads_${FILE}_percent_bad.fasta | python ./evaluate.py
 echo
 echo
done
