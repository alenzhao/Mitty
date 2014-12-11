Mitty is a collection of programs that enable us to generate simulated genomic data to test our algorithms.

Mitty generates simulated data by starting from a reference genome and then adding controlled mutations to it to create
a variant "sample" genome. Mitty can also simulate populations of genomes derived via sexual reproduction. Mitty can
then generate reads from the chromosomes of the mutated genome(s). 

At each stage of this process we know the origins of the data and so we know

- The correct alignment of each read
- The correct variants in the "sample" that these reads came from.

**The data generated by Mitty can, therefore, be used to test the correctness of aligners and variant callers under
different conditions and data set characteristics.**

Mitty, at its heart, is a set of command line programs that are called with various parameters to generate different 
simulated data. Mitty has, additionally, been wrapped for the Seven Bridges Genomics computation platform (IGOR) and
can be used from that platform to generate simulated data.

The important components of Mitty are described as::

    Reference Genome  + denovo.py      -->  .vcf file
        mutations of various kinds representing a single sample
  

    Reference Genome  + population.py  -->  lot of .vcf files
        a founder population of genomes
      

    Reference Genome + founder population  + generations.py  -->  lot of .vcf files, lineage files
        multiple generations created by sexual reproduction
      
      
    .vcf file + vcf2reads.py           -->   .fastq file
        perfect and corrupted reads from "sample"
    
    .fastq file + reads2bam.py         -->   .bam file
        perfect alignment of "sample" reads to original reference
    
    
    .bam file + checkbam.py            -->  diagnostic files
        Allows us to check the quality of our alignments

Installation
============

Mitty can be installed using pip like this::

    pip install git+ssh://git@gitlab.sbgenomics.com:graphgenome/mitty.git#egg=mitty


This requires that you are connected to VPN, and that you have an SSH key setup on Gitlab.


Testing
=======
You should have ``pip install nose`` previously. Then simply run ::

    nosetests mitty


Developing
==========

To develop on Mitty, simply clone the repository, and from the project root run ::

    pip install -e .
