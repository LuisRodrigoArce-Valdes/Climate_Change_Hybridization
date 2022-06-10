#!/bin/sh
# 07_COIs_Meta.sh
# By Luis Rodrigo Arce Valdés (09/06/22)
# Using this script we will create a database of the meta information of all COI sequences used in this project

# Creating files with the meta information of both BOLD and GenBank sequences
# BOLD
mkdir -p ../data/08_Meta
grep "COI-5P" ../data/04_Muscle/* > ../data/08_Meta/Bold.txt
# GeneBank
grep "|gb|" ../data/04_Muscle/* > ../data/08_Meta/GB.txt



