#!/bin/bash

# set path
GENOM = CDC1
DIR = /home/VBarton/Documents/University/Smajs/SeqData/$GENOM

pools="POOL1 POOL2 POOL3 POOL4 TprC TprE"

################
cd $DIR
mkdir vcf

for i in $pools
do
  # dir to poolX
  cd $i

  # assembly
  bwa index ref.fa
  bwa mem -MT 30 -t 16 ref.fa Reads1.fastq Reads2.fastq > aln.sam

  # cleaning
  samtools faidx ref.fa
  samtools view -bT ref.fa -@ 15 -o pool.bam aln.sam
  samtools sort -o pool_srt.bam -@ 15 pool.bam
  samtools index pool_srt.bam

  samtools rmdup pool_srt.bam pool_rm.bam
  samtools index pool_rm.bam

  samtools view -b -F 4 -@ 15 -o pool_clear.bam pool_rm.bam
  samtools index pool_clear.bam

  vcf = "${i}_stat.vcf"
  samtools mpileup -f ref.fa -Q 50 -d 200000 -uvI -t AD,DP,INFO/AD -o $vcf pool_clear.bam

  cp $vcf ../vcf/$vcf
done
