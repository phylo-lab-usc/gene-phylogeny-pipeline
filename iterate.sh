#! /usr/bin/bash

echo "type in gene"
read geneID
touch $geneID.fasta

FILES="calbicans.fasta
anidulans.fasta
afumigatus.fasta
scerevisiae.fasta"

for f in $FILES
do
      awk -v pat="$geneID" '$0~pat' $f >> $geneID.fasta 
done
