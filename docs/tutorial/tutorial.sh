  genomes
  genomes --version
 genomes show parameters
 genomes show model-list
  genomes show variant-model snp
  genomes show spectrum-model double_exp
 genomes show population-model standard
 genomes generate --dry-run variations.json
 genomes generate variations.json -v -p
 genomes genome-file summary reddus_genomes.h5
 genomes genome-file sfs reddus_genomes.h5 1
 genomes genome-file write-vcf reddus_genomes.h5 g0_s0.vcf --sample-name 'g0_s0'
 head -n 20 g0_s0.vcf
 reads show parameters
 reads show model-list
 reads show read-model simple_illumina
  reads generate reads.json -v
  reads --qname
 head -n 16 reads_c.fq | awk '/^@/ {print}'
  cat reads.bed
 bwa index reddus_pentalgus.fa.gz
 samtools tview -d T -p "NC_010142.1:50" bwa.bam reddus_pentalgus.fa
 perfectbam --help
 perfectbam --tags
  perfectbam -v -v bwa.bam
  perfectbam bwa_poor.bam
 misplot --help
 misplot bwa_bad.bam --matrix bwa_mat.png --circle bwa_cir.png
 misplot bwa_poor_bad.bam --matrix bwa_poor_mat.png --circle bwa_poor_cir.png
 indelbam --help
 indelbam --tags
 indelbam bwa_per.bam reddus_genomes.h5 g0_s0 bwa_indel.bam bwa_indel.pkl bwa_indel_summary.json --indel-range 50
 indelbam bwa_poor_per.bam reddus_genomes.h5 g0_s0 bwa_poor_indel.bam bwa_poor_indel.pkl bwa_poor_indel_summary.json --indel-range 50
 cat bwa_indel_summary.json
 indelplot --help
 indelplot --title 'BWA' bwa_indel.pkl bwa_indel.png --indel-range 100
 mqplot --help
 mqplot bwa_per.bam --title 'BWA' -o bwa_mq.png
 aligner-summary-report --help
 aligner-summary-report bwa g0_s0 bwa_indel_summary.json novel-bwa_indel.png known-bwa_indel.png bwa_mq.png bwa_cir.png bwa_mat.png report.html --db-summary reddus_genomes_summary.txt
 badbams --help
 badbams bwa_bad.bam bwa_poor_bad.bam -v -p
 acubam --help
.. #  command-output:: acubam bwa.bam reddus_genomes.h5 --sample-name "g0_s0" -p
 genomes from-vcf --help
 genomes from-vcf g0_s0.vcf.gz converted.h5 --sample-name g0_s0 -v
 genomes genome-file summary converted.h5
