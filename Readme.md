Mitty is a collection of modules and scripts that enable us to generate simulated genomic data to test our algorithms.


User guide
==========
There are two branches in the repository:
    **master** - stable working code
    **dev**    - code could be unstable/unworking but will have the latest experimental stuff going on

The different components of the simulator are
mutate.py -

Example workflow
----------------

## Test the read generator with ideal reads from a reference genome
```
python ~/Code/Mitty/reads.py --ref=porcine_circovirus.fa --read_len=50 --read_count=100  #Generate the reads
bwa index porcine_circovirus.fa
bwa mem porcine_circovirus.fa simulated_reads.fastq > aligned.sam
samtools view -Sb aligned.sam > aligned.bam
samtools sort aligned.bam aligned_sorted
samtools index aligned_sorted.bam
samtools tview aligned_sorted.bam porcine_circovirus.fa
```

## Test the SNP generator with ideal reads from a mutated genome
```
python ~/Code/Mitty/mutate.py --ref=porcine_circovirus.fa --out=mutated --paramfile=params.py  #Generate SNPs
python ~/Code/Mitty/reads.py --ref=mutated.fa --read_len=50 --read_count=100  #Generate the reads
bwa index porcine_circovirus.fa
bwa mem porcine_circovirus.fa simulated_reads.fastq > aligned.sam
samtools view -Sb aligned.sam > aligned.bam
samtools sort aligned.bam aligned_sorted
samtools index aligned_sorted.bam
samtools tview aligned_sorted.bam porcine_circovirus.fa   # You should see the SNPs marked out
samtools mpileup -uf porcine_circovirus.fa aligned_sorted.bam | bcftools view -bvcg - > var.raw.bcf
bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 > var.flt.vcf
# You should be able to compare var.filt.vcf and mutated_variants.vcf and verify they are the same
```

Where `params.py` is the same as the example parameter file in the repository

```
"""Example parameter file for mutate program

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
snp = {
  'p': 0.01
}
```


Dev notes
=========
### Parameter files for mutate
1. I chose to use parameter files because we often want to rerun experiments and it became clear early on that there would
be a lot of parameters.
1. I chose to use python for the parameter file for parsimony and flexibility
1. The parameter distribution between file and command line was based on predictions of which parameters we could
experiment with most during testing

(Since you were dying to know: Mitty comes from James Thurber's "The Secret Life of Walter Mitty" one of my favourite
pieces from one of my favourite authors. Though [Wikipedia][wiki] has a less favourable interpretation of what Walter Mitty
stands for I follow the interpretation found in the 2013 movie of the same name. Life is difficult and full of
insurmountable obstacles. If you do not even dream that you have surmounted these obstacles how are you going to even
start?)

[wiki]: http://en.wikipedia.org/wiki/The_Secret_Life_of_Walter_Mitty