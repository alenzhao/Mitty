Tutorial
========

We will introduce Mitty via a data simulation problem. Let us say our colleague Heng Li has written an aligner
called ``bwa`` and a variant caller called ``samtools`` and we want to figure out how well they do using synthetic data.

Our plan is as follows:

* Consider a small reference genome (we will use a made up genome *Reddus pentalgus* made up of the five largest
chromosomes of the Red Alga, *Cyanidioschyzon merolae strain 10D*, record_, ftp_ and pretend it's a diploid organism)
* Sprinkle some variants on the reference genome to create a population of 1000 sample genomes (why not?)
* Take 30x Illumina like reads from one of the sample genomes
* Use `bwa` to align the reads back to the reference genome
* Assess how well `bwa` aligned the reads


.. _record: http://www.ncbi.nlm.nih.gov/assembly/GCF_000091205.1/#/def
.. _ftp: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000091205.1_ASM9120v1

--------------

For the impatient, the commands we'll be using for this experiment are:

.. literalinclude:: tutorial/tutorial.sh

This file is included in the source distribution under ``docs/tutorial/tutorial.sh``

--------------

Creating genomes
----------------
The example reference genome has 20 chromosomes. We'd like to add SNPs, insertions and deletions to chromosome 1 and 20.
Since this is for testing purposes we would like to have a large number of indels distributed across a range of lengths,
rather than the power-law distribution we see in nature (which gives us many short indels and far fewer longer indels).

To perform this simulation, we need to create a parameter file describing our simulation requirements and pass this to
the ``genomes`` program. In the parameter file we need to describe three components:

 * the individual variant models
 * the site frequency spectrum model and,
 * the population generation model

The variant models are used to sprinkle variations (mutations) on top of an input reference genome. These variations are
put into a master list of variations. Each stock variation model has a parameter, ``p``, that sets the
per base probability of observing the variation being generated. For example, setting ``p=0.001`` for the SNP plugin
and ``p=0.0001`` for the delete plugin will give us SNP and delete densities that are similar to that observed in typical
human samples.

Internally, the plugins scale the probability according to the site frequency spectrum such that samples
drawn from the master list contain variations obeying the site frequency spectrum and, overall, the desired ``p`` value.

The population generation model specifies how to construct individual samples from the master list of variants. The
default population model randomly selects variations for each sample from the master list following the site frequency
spectrum. This default model is sufficient for large population simulations.

For greater detail please see [theory and algorithms] XXX

The genomes program performs several tasks. Typing ``genomes`` on the command line will produce a list of sub-commands.

.. command-output::  genomes

For the task at hand we need the ``genomes generate`` tool. We need to construct a parameter file which contains,
among other things, specifications for variant models, spectrum models and population models.

.. command-output:: genomes show parameters

What are the models available to us?

.. command-output:: genomes show model-list

What kind of parameters does, for example, the SNP model need?

.. command-output::  genomes show variant-model snp

The spectrum model?

.. command-output::  genomes show spectrum-model double_exp

The standard population model?

.. command-output:: genomes show population-model standard

Let us create a parameter file for the simulations, calling it ``variations.json``:

.. literalinclude:: tutorial/variations.json
    :language: json
    :linenos:

*Note that if we are happy with a default plugin parameter, we can omit that in the parameter file.*

Before we run the actual simulation, let's make sure the site frequency spectrum looks satisfactory

.. command-output:: genomes generate --dry-run tutorial/variations.json

Let's run this command and create a database of simulated genomes

.. command-output:: genomes generate tutorial/variations.json -v

(*Using the -v option will give detailed logger messages, using the -p option will give a progress bar*).

Let's take a peek at the produced genome population database:

.. command-output:: genomes genome-file summary tutorial/red_alga_genomes.h5

We can also take a look at the site frequency spectrum in the generated population, for chromosome 1, for example

.. command-output:: genomes genome-file sfs tutorial/red_alga_genomes.h5 1

Right now, all our genomes are in a database. The rest of the world works in VCF files, so let's write out one of the
samples as a VCF file

.. command-output:: genomes genome-file write-vcf tutorial/red_alga_genomes.h5 tutorial/red_alga --sample-name 'g0_s0'
.. command-output:: head -n 30 tutorial/red_alga.vcf | tail -n 10
    :shell:

Things seem to have run satisfactorily. Note that we produce phased VCF files and we also use the notation `1|0` since we
have complete knowledge of the phasing. Now let's generate a bag of reads from this genome.

Taking reads
------------
We will use the `reads` program to generate a fastq file of reads from the selected genome. We could do with a refresher
on how to create a parameters file for `reads`:

.. command-output:: reads show parameters

And what sort of read models are currently available?

.. command-output:: reads show model-list

``simple_illumina`` sounds promising, what kind of parameter file does it need?

.. command-output:: reads show read-model simple_illumina

Ok, let's put together a parameter file for our experiment, let's call it `reads.json`:

.. literalinclude:: tutorial/reads.json
    :language: json
    :linenos:

Let's run this command and create some reads from the variant genome

.. command-output::  reads generate tutorial/reads.json -v

Testing alignment accuracy of BWA MEM
-------------------------------------
First, use BWA-MEM to create an alignment:

.. command-output:: bwa index tutorial/red_alga.fa.gz
    :ellipsis: 1,-2
.. command-output:: bwa mem -t 8 -p tutorial/red_alga.fa.gz  tutorial/reads_c.fq > tutorial/temp.sam
    :shell:
    :ellipsis: 1,-2
.. command-output:: samtools view -Sb  tutorial/temp.sam > tutorial/temp.bam
    :shell:

.. command-output:: samtools sort tutorial/temp.bam  tutorial/reads
.. command-output:: samtools index tutorial/reads.bam

Let's check the alignments:

.. command-output:: samtools tview -p NC_010127.1:350 -d T tutorial/reads.bam tutorial/red_alga.fa.gz


We can use an alignment browser, such as [Tablet]_, to see the alignments

.. image:: _static/reads_bam_tablet.png

(You can see that the coverage is what we dialled in. The one red colored base is a read error - since we are inspecting
the corrupted reads file)

.. [Tablet] Milne I, Stephen G, Bayer M, Cock PJA, Pritchard L, Cardle L, Shaw PD and Marshall D. 2013. Using Tablet for visual exploration of second-generation sequencing data. Briefings in Bioinformatics 14(2), 193-202. (`web`_)
.. _web: http://ics.hutton.ac.uk/tablet/


Analyse the alignment
.....................

We use ``perfectbam`` to both process_bams the performance of the aligner and to create a perfectly aligned ``.bam`` file from the
``.bam`` file output by the aligner.

.. command-output::  perfectbam ../examples/demo/Out/reads.bam


`perfectbam` creates a `.json` file with some summary performance statistics and a `.csv` file with details about the
misaligned reads.

.. command-output:: head -n 15 ../examples/demo/Out/reads_summary.json

.. command-output:: head -n 10 ../examples/demo/Out/reads_misaligned.csv

We all know how exciting it is to scroll through a .csv file, so we have some additional tools to analyse the misaligned
reads and display the information graphically (These require Matplotlib_ to be installed).

*(For aesthetic reasons, we will use a different reads file that has reads from all the chromosomes, but at lower coverage
so it doesn't take so long to run)*

.. command-output:: ../examples/demo/low_coverage_reads.sh
  :ellipsis: 0

.. command-output:: plot_align circle ../examples/demo/Out/null_reads_low  --down-sample 2
.. image:: ../examples/demo/Out/null_reads_low_circle_plot.png

This is a plot of where the reads end up after alignment. As can be seen, while a small percentage of reads end up in
different chromosomes, most of the misalignments are local (the lines end up on the same chromosome).

.. command-output:: plot_align matrix ../examples/demo/Out/null_reads_low
.. image:: ../examples/demo/Out/null_reads_low_matrix_plot.png


.. _Matplotlib: http://matplotlib.org/index.html

Taking reads just from regions around variants
----------------------------------------------
By setting the ``variants_only`` flag to *True* in the parameters file we can selectively generate reads from regions bordering
a variation.

.. AArgh. Could not get Sphinx to highlight the line (https://github.com/ClusterHQ/flocker/issues/337) - I think the lame green highlight matches the background

.. literalinclude:: ../examples/demo/illumina_reads_variants_only.json
    :language: json
    :emphasize-lines: 15
    :lineno-match:
    :lines: 13-17

.. command-output:: reads generate  ../examples/demo/illumina_reads_variants_only.json -v
.. command-output:: bwa mem -t 8 -p  ../examples/data/red_alga.fa.gz  ../examples/demo/Out/vreads_c.fq > ../examples/demo/Out/temp.sam
    :shell:
    :ellipsis: 1,-2
.. command-output::  samtools view -Sb  ../examples/demo/Out/temp.sam > ../examples/demo/Out/temp.bam
    :shell:

.. command-output:: samtools sort  ../examples/demo/Out/temp.bam  ../examples/demo/Out/vreads
.. command-output:: samtools index  ../examples/demo/Out/vreads.bam

We can use an alignment browser, such as Tablet, to see the alignments

.. image:: _static/vreads_bam_tablet.png

Note how the reads are restricted to the variant locations, as we asked for.


Analysing variant callers
-------------------------
TODO