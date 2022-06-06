#!/bin/sh
# 01_Getting_COIs.sh
# By Luis Rodrigo Arce Valdés (05/06/22)
# Using this script we will download COI sequences from the BOLD database hybridizing pair species
# http://www.boldsystems.org/
# https://github.com/CNuge/BOLD-CLI

# 01.- Copying database
mkdir -p ../data/01_raw/
cp ../../01_Review/results/Database.txt ../data/01_raw/

# Changing apis mellifera "subspecies" for the italian one
sed -i 's/Apis mellifera ssp/Apis mellifera ligustica/g' ../data/01_raw/Database.txt
echo "We used Apis mellifera ligustica as representative of one of the subspecies that have hybridized with the africanized bees" > ../results/Bees_README.md

# 03. Creating output diretories
mkdir -p ../data/02_BOLD

# 04.- First I am doing a list of all hybridizing species
# Second column (First species):
cut -f2 ../data/01_raw/Database.txt | tail -n +2 > tmp1

# Third column (Second species):
cut -f3 ../data/01_raw/Database.txt | tail -n +2 > tmp2

# Concatenating both columns, sorting, removing spaces at the end of strings and removing duplicates
cat tmp1 tmp2 | sort | uniq | sed 's/ *$//' > ../data/02_BOLD/Species.txt
rm tmp1 tmp2

# 05.- Searching COI sequences within the BOLD systems database
./bold-cli -output ../data/02_BOLD/COIs.fasta -query sequence -marker COI-5P -taxon ../data/02_BOLD/Species.txt
