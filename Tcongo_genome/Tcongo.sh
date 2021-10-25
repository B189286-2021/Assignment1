#!/bin/bash

#Create a Tcongo_genome directory in Assignment1 folder
mkdir ~/Assignment1/Tcongo_genome

#Copy in the genome sequence
cp /localdisk/data/BPSM/AY21/Tcongo_genome/*Genome.fasta.gz ~/Assignmet1/Tcongo_genome

#Unzip the file
gunzip *gz

#Create an index for the Trypanosoma brucei brucei reference genome
bowtie2-build ~Assignment1/Tcongo_genome/*Genome.fasta Tcongo_mapping.btindex

#review folder
ls

#Change directory to the fastq folder - fastq files to be read
cd ~/Assignment1/fastq

#Riview folder
ls

#Copy in the indexed files
cp ~/Assignment1/Tcongo_genome/*btindex* ~/Assignment1/fastq

#Map them reads - output is sam file then convert to bam with samtools view
for i in *1.fq
do
base=$(basename $i "_1.fq")
bowtie2 -p 4 -x ~/Assignment/fastq/Tcongo_mapping.btindex -1 ${base}_1.fq -2 ${base}_2.fq | samtools view -b -o {base}.bam
done

#sort bam files using samtools sort
for i in *bam 
do
base=$(basename $i ".bam")
samtools sort -O bam -o ${base}.sorted.bam ${base}.bam
done

#index files using samtools index
for i in *sorted.bam
do
base=$(basename $i ".sorted.bam")
samtools index {base}.sorted.bam
done

#Create a directory and copy in the "bedfile" that contains the information about the gene location
mkdir ~/Assignment1/bedtools
cp /localdisk/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed

#generate counts using bedtools multicov
cd ~/Assignment1/fastq/
for i in *.sorted.bam
do
base=$(basename $i ".sorted.bam")
bedtools multicov -s -bams ${base}.sorted.bam -bed ~/Assignment1/bedtools.TriTrypDB-46_TcongolenseIL3000_2019.bed > ${base}_counts.bed
done


#!usr/bin/awk -f 
#count no. of overlapping features
for i in *counts.bed
do
base=$(basename $i ".counts.bed")
cat ${base}counts.bed |awk '{
FS="\t" ;
if(NF == 6)
{print $0 ;}
}' | wc -l
done

#count no. of feature with non-zero counts
for i in *counts.bed 
do
base=$({base}.counts.bed)
cat {base}.counts.bed |\
> awk '{if($6> 0) print $6'} | wc -l
done


