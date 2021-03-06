#!/usr/bin/env bash
set -x


mitty -v4 bam2illumina sample.bam ./rd-model.pkl "This model is taken from a BAM of unknown provenance" --every 5 --min-mq 30
mitty describe-read-model ./rd-model.pkl model.png


mitty -v4 filter-variants hg001.vcf.gz INTEGRATION hg001.bed - 2> filter.log | bgzip -c > hg001.filt.vcf.gz
tabix -p vcf hg001.filt.vcf.gz


mkfifo tf1
mkfifo tf2
mitty -v4 generate-reads ~/Data/human_g1k_v37_decoy.fasta hg001.filt.vcf.gz INTEGRATION hg001.bed ./rd-model.pkl 30 7 >(tee > tf1 >(gzip > r1.fq.gz)) --fastq2 >(tee > tf2 >(gzip > r2.fq.gz)) --threads 2 &
mitty -v4 corrupt-reads ./rd-model.pkl <(cat < tf1) >(gzip > r1c.fq.gz) 7 --fastq2-in <(cat < tf2) --fastq2-out >(gzip > r2c.fq.gz)

bwa mem ~/Data/human_g1k_v37_decoy.fasta r1c.fq.gz r2c.fq.gz | samtools view -bSho out.bam
samtools sort out.bam > bwac.bam
samtools index bwac.bam

mitty -v4 god-aligner ~/Data/human_g1k_v37_decoy.fasta r1c.fq.gz perfect.bam --fastq2 r2c.fq.gz --threads 4

bam diff --in1 bwac.bam --in2 perfect.bam --out diffd.unsorted.bam
samtools sort diffd.unsorted.bam > diffd.bam
samtools index diffd.bam

# Cleanup
# rm *.bam *.bai *.png
