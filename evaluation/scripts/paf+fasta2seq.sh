#!/bin/bash

PAF=$1
FASTA=$2

cat $PAF | awk -v OFS='\t' '{print $1, $3, $4, "", "", $5}' > $PAF.query.bed
cat $PAF | awk -v OFS='\t' '{print $6, $8, $9, "", "", "+"}' > $PAF.target.bed

bedtools getfasta -fi $FASTA -bed $PAF.query.bed -s | grep '^>' -v | sed 's/^/>/' > $PAF.query.tmp
bedtools getfasta -fi $FASTA -bed $PAF.target.bed -s | grep '^>' -v | sed 's/^/</'> $PAF.target.tmp

paste -d '\n' $PAF.query.tmp $PAF.target.tmp > $PAF.seq

rm $PAF.query.bed $PAF.target.bed $PAF.query.tmp $PAF.target.tmp
