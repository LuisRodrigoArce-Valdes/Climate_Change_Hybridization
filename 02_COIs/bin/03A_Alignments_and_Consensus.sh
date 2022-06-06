#!/bin/sh
# 03_Alignments_and_Consensus.sh
# By Luis Rodrigo Arce Valdés (28/03/22)
# With this script we will align all COI sequences for each species and then create a Consensus sequence for each one.
# To use for de novo analyses, without saved alignemts and consensus sequences
mkdir -p ../data/04_Muscle
mkdir -p ../data/05_Consensus
for n in $(ls ../data/03_Species_Fastas/)
do
	echo "########"
	echo "$n"
	grep ">" ../data/03_Species_Fastas/$n | wc -l
	# Aligning with muscle's "super5" algorithm (faster than the original method in large databases)
	muscle -super5 ../data/03_Species_Fastas/$n -output ../data/04_Muscle/$n.afa
	# Creating consesus sequence [cons only works if we have more than two sequencues]
	# Using -plurality = 1, with at least one sequence per position consensus sequence will have a called nucleotide
	if [ $(grep ">" ../data/03_Species_Fastas/$n | wc -l) -gt 1 ]
	then
		cons -sequence ../data/04_Muscle/$n.afa -outseq ../data/05_Consensus/$n -name $n -plurality 1 -verbose
	else
		cp ../data/04_Muscle/$n.afa ../data/05_Consensus/$n
		sp=">$n"
		sed -i "1s/.*/$sp/" ../data/05_Consensus/$n
	fi
done
