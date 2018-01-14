for FILE in 1 1_to_3 3_to_5
do
 echo ${FILE}
  python assembly.py ./sample_data/reads_${FILE}_percent_bad.fasta ./sample_data/output_${FILE}_percent_bad.fasta
 echo
 echo
done


