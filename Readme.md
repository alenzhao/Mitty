Mitty is a collection of modules and scripts that enable us to generate simulated genomic data to test our algorithms.
The scripts allow us to simulate mutations on a reference sequence/genome and then simulate reads from that mutated
sequence/genome.

Quickstart
----------

The process for running Mitty components to create reads from a mutated genome starting from only a reference
sequence is illustrated schematically below. Please see `complete_example.sh` under the `Recipes` directory
(and the relevant parameter files under the `Params` directory) to understand the details of the command line
invocations and parameter files. For advanced users please see the rest of the docs and the `Plugins` directory
to understand how to write python code to simulate different kinds of mutations and reads.

                    ----------
         fasta     |          |---> smalla file
          file --->| converta |
                   |          |---> heada file
                    ----------

For efficiency purposes we strip the original fasta file of the header and all new lines. The
       resulting file is called a smalla file and is what the rest of the tools use. This conversion
       can be done easily using the converta.py script. The header is saved into a .smalla.heada file


                    mutation
                   parameters
                       |
                       V
                    --------
                   |        |
       ref seq --->| mutate |----> VCF file
                   |        |----> side car file with sim params
                    --------

Given a set of mutation instructions we can use mutate.py to generate a VCF file. For further
       processing the vcf file should be compressed using `bgzip` and indexed by `tabix`.



                    ---------
       ref seq --->|         |
                   | vcf2seq | ---> mutated seq
       VCF     --->|         |
                    ---------

Using the vcf2seq tool we can write out the mutations indicated by VCF into a complete mutated
       sequence saved as a smalla file.


                     read
                   parameters
                       |
                       V
                    --------
                   |        |----> corrupted ("real") reads (BAM/FASTQ)
                   |        |
           seq --->| reads  |----> ideal reads (BAM/FASTQ)
                   |        |
                   |        |----> side car file with sim params
                    --------

The reads tool enables us to take a sequence and generate simulated reads from it. The reads can
       simulate various error and property profiles of different sequencers





Each module is designed to run as a script. Typing `python mutate.py -h` or simply `python mutate.py` etc. will list
usage and input requirements.

There are two branches in the repository:

    master - stable working code
    dev    - code could be unstable/unworking but will have the latest experimental stuff going on

The code requires the following non-standard modules

    BioPython   - pip install biopython --user
    PyVCF       - pip install pyvcf --user
    pysam       - pip install pysam -- user

The code requires the following external tools to run the examples

    bgzip       - to compress VCF file
    tabix       - to index VCF file
    samtools    - indexing fasta and converting between bam and sam etc
    bwa         - alignments etc.



Subdirectories
--------------
    Params      - example parameter files for mutate.py and reads.py
    Recipes     - snippets of code (shell scripts and python scripts) to do/show particular tasks. useful for devs and
                  users alike
    Data        - test data for the programs
    Plugins     - directory where simulation models are stored

Data
----
 - porcine_circovirus.fa - 702bp (http://www.ncbi.nlm.nih.gov/nuccore/AY735451.1)
 - adenovirus.fa   -  34094bp  (http://www.ncbi.nlm.nih.gov/nuccore/AB026117.1)


Files and formats
-----------------

## VCF

`vcf2seq.py` currently handles a specific interpretation of the VCF. In the most liberal interpretation of the VCF the
REF sequence corresponds to bases matching the reference starting at POS and those bases are replaced by the bases found
in ALT. Any number of bases in ALT may match the REF in any place (though this is not very useful to us).

Mitty has a more strict interpretation of the VCF. In Mitty's interpretation, at most, the first base
of REF will match with ALT. All other bases must be different. All variants can be coded in this manner:

Say our original sequence is `ATCGGATC`

                 POS   REF     ALT
    Deletion      1   ATCG     A     -> AGATC
    Deletion      2   TCG      .     -> AGATC     (Though this form is interpreted by vcf2seq it is never produced by mutate.py)
    Deletion      1    A       .     -> TCGGATC
    Insertion     0    .      TTT    -> TTTATCGGATC
    Insertion     1    A      AGGG   -> AGGGTCGGATC
    SNP           1    A       G     -> GTCGATC

## "Buffer bases" between simulated variants

If we have variants adjacent to each other the most parsimonious description of the resulting variation can be different
from the original variants.

For example, considering `M=ATCGATCG` and an insertion and deletion as follows

    POS REF ALT
    1   A   ACC
    1   AT  A

We get

         1  2345678
    R    A  TCGATCG
    M    ACC CGATCG

This can, actually, be most parsimoniously expressed as

         1 2345678
    R    A TCGATCG
    M    ACCCGATCG

Which is a single base insertion followed by a SNP

    POS REF ALT
    1   A  AC
    2   T  C

For reasons of such ambiguity Mitty places a minimum 1 base "buffer" between variants, making the generated variant
identical to the most parsimonious description.


Dev notes
=========

Plugin system for simulation models
-----------------------------------

Mitty implements models for variant and read simulation as Python modules located in the Plugins directory. The modules
need to expose a few key functions that `mutate.py` and `reads.py` use to determine variant and read characteristics.
Mitty infers the module name from the name given in the parameter .json file. Mitty comes with some stock models for
variant and read simulation which can be used as example code.

### Random number generators, blocked computation etc.
You will note that the stock plugins accept one or more inputs that serve as seeds for internal random number
generators. These seeds need to be specified in the simulation parameter files and ensure reproducibility of
simulations.

Some of the models use several, independently seeded, random number generators for different variables (e.g.
insertion length and insertion position) to avoid unexpected interactions between such simulated variables. The
algorithms are also designed such that the simulation results are independent of block size.


Algorithms
----------
### Variants, reads and CIGARS
One big goal of Mitty is to serve up realistic test data for bioinformatics algorithms, from aligners to variant
callers. Testing whether a variant caller is correctly working on the simulated data is relatively easy: we simply
compare the variant caller's VCF file with the answer book VCF generated by Mitty. It is, however, slightly more involved
to deduce if an aligner is correctly aligning the simulated reads, and to diagnose how the performance of an aligner is
affecting the accuracy of a variant caller. To this end Mitty has a system to compute the correct read position and
CIGAR for each simulated read. This information is stored in the read's qname string so that it is easily accessible
to diagnostic programs.

In order to generate reads based on a given VCF file and a reference sequence we go through a two step process.
We first generate the mutated sequence (`mut_seq`) along with some other information that encodes the difference between
each base in the `mut_seq` and corresponding positions on the `ref_seq`. We then generate reads from the `mut_seq` using
the sidecar information to compute the correct POS values and CIGAR strings for the reads.

The algorithm is best introduced through a series of examples. In the examples the reference sequence is labelled `R` and
the mutated sequence is labeled `M`. The information for setting the POS and CIGAR for the read is taken from an
array `pos` that accompanies `M`

#### Generating `pos`

Let `R = ACTGACTG`

Consider a single base insertion at position 1

    POS REF ALT
    1   A   AT

         1 2345678
    R    A CTGACTG
    M    ATCTGACTG
    pos  1223456789


Consider a multiple base insertion at position 1

    POS REF ALT
    1   A   ATT

         1  2345678
    R    A  CTGACTG
    M    ATTCTGACTG
    pos  12223456789


Consider a multiple base insertion at last position

    POS REF ALT
    8   G   GTT

         12345678
    R    ACTGACTG
    M    ACTGACTGTT
    pos  12345678999

Consider a multiple base deletion

    POS REF ALT
    2   CTG  C

         12345678
    R    ACTGACTG
    M    AC  ACTG
    pos  12  56789

Consider a SNP, an insertion and a deletion

    POS REF ALT
    2   C   T
    4   G   GTT
    6   CTG C

         1234  5678
    R    ACTG  ACTG
    M    ATTGTTAC
    pos  123455569


`pos` is generated by copying over the index from `R`. When we encounter an insertion we copy over the index of the next
reference base as many times as there is an insertion. Deletions are simply skipped. For the purposes of computing `pos`
we also add an imaginary base position at the end of the reference sequence (9 in this case)

#### Generating CIGARS and POS for reads from `pos`
Consider our last example and some reads from `M`

         1234  5678
    R    ACTG  ACTG
    M    ATTGTTAC
    pos  123455569
         ++++---------> POS = 1 (The first pos value we encounter)
                        CIGAR = 4M  (2-1=1 -> 1M
                                     3-2=1 -> 1M
                                     4-3=1 -> 1M
                                     5-4=1 -> 1M)

    M    ATTGTTAC
    pos  123455569
          ++++--------> POS = 2 (The first pos value we encounter)
                        CIGAR = 3M1I  (3-2=1 -> 1M
                                       4-3=1 -> 1M
                                       5-4=1 -> 1M
                                       5-5=0 -> 1I)

    M    ATTGTTAC
    pos  123455569
           ++++-------> POS = 3
                        CIGAR = 2M2I  (4-3=1 -> 1M
                                       5-4=1 -> 1M
                                       5-5=0 -> 1I
                                       5-5=0 -> 1I)

    M    ATTGTTAC
    pos  123455569
             ++++-----> POS = 5
                        CIGAR = 2I2M  (5-5=0 -> 1I
                                       5-5=0 -> 1I
                                       6-5=1 -> 1M
                                       9-6=3 -> 1M + 2D) The D only comes into play if our read crosses the deletion

To see how a deletion affects our POS and CIGAR consider another previous example

    POS REF ALT
    2   CTG  C

         12345678
    R    ACTGACTG
    M    AC  ACTG
    pos  12  56789
         ++  ++-------> POS = 1
                        CIGAR = 2M2D2M  (2-1=0 -> 1I
                                         5-2=3 -> 1I + 2D The 2D comes into play because the read crosses the boundary
                                         6-5=1 -> 1M
                                         7-6=1 -> 1M)

Example of an unmapped read

    POS REF ALT
    2   C  CAATTGG

         12      345678
    R    AC      TGACTG
    M    ACAATTGGTGACTG
    pos  123333333456789
           ++++-------> POS = 3
                        CIGAR = 4I  (3-3=0 -> 1I
                                     3-3=0 -> 1I
                                     3-3=0 -> 1I
                                     3-3=0 -> 1I)
    For a read to be mapped, there has to be at least one M. Since there are no Ms we discard the POS and CIGAR as this
    is an unmapped read

`reads.py` generates simulated reads from `mut_seq` based on the read model. Using the `pos` arrays it
also generates appropriate alignment information (POS and CIGAR) that is stored in the qname string.
(Note that while the BAM specs do not place a limit on the length of the qname string both Tablet and IGV expect a
string with length < 255 characters. It is possible that the qname will exceed this and you won't be able to open a
set of simulated reads using tools that arbitrarily limit the qname). If no `pos` file is supplied `reads.py` assumes
we are taking reads from a reference sequence and the POS values are actual positions of the reads and all the cigars
are of the form `100M` (For e.g. 100 base reads).

Computing POS: For every read, the POS value is simply the index from `pos` corresponding to the first base of the read
EXCEPT for unmapped reads.

Computing the CIGAR:

1. Initialize the base counter to `None`, set mapped flag to `False`
2. Step through the each base of the read and look at the difference in `pos` values `dp`
3. If `dp==1`, if the counter is any thing other than `M`, flush it. Set or increment counter as `M`. Set mapped flag to `True`
4. If `dp==0`, if the counter is other than `I`, flush it. Set or increment counter as `I`
5. If `dp>1`, if the counter is other than `M`, flush it. Set and flush counter as `M`, set counter as `D` to be dp-1
6. Continue from 2 until done.
7. Flush any counter other than `D`
8. If the mapped flag is `False` reset POS and CIGAR - this is an unmapped read.


Misc design choices
--------------
### Choice to output a mutated sequence as a whole
Though generating a whole mutated sequence uses a lot of disk space, I chose this approach as it ended up being simpler
than coming up with an algorithm for generating reads on the fly based ona  VCF file. In the future the code may be
converted to do on the fly generation.

### Parameter files for mutate
1. I chose to use parameter files because we often want to rerun experiments and it became clear early on that there would
be a lot of parameters.
1. I chose to use python for the parameter file for parsimony and flexibility
1. The parameter distribution between file and command line was based on predictions of which parameters we could
experiment with most during testing
1. At this time I do not know whether having everything on the commandline would be better for PIPITOR or if param files
are preferred for Platform integration, but either way is a short code reorganization that can be done quickly at the
time of integration.

### POS and DIFFPOS files
These are simple binary files carrying unsigned 4 byte int information. This is enough to handle index/index diff sizes
for human genome sizes, though if we ever work on heavily mutated specimens of the loblolly pine, perhaps we have to
go to 8 byte ints.



Disorganized - don't read below this
====================================



### Generating reads from a ref seq and VCF(s)
We want to be mindful that genomic sequences can be quite large and we might not want to carry multiple copies of such
sequences around, especially when we are simulating reads from heterogenous sequences. For this reason we choose to
avoid generating the complete sequences or ever loading the whole ref seq into memory at once.

The conceptual way to generate reads is as follows:

1. Generate the variant sequences as a whole -> there may be more than one variant sequence depending on whether there
are multiple variants at the same locus. The number of variant sequences depends on the combinatorials of all the
multiple variants
2. Take sections from each sequence at random locations, according to a predetermined distribution

If we want to avoid physically generating all combinations of variant sequences we do the following.

1. We move along the ref-seq in blocks
2. We consider only the variants within that block
3. We generate all variant sequences for that block and then generate reads randomly within that block

Notes:
1. How to handle non-local mutations? (translocation etc.)
2. How to handle multiple variants - how to link variants in one mutation with variants in another? - will currently do
a simple combinatorial with equal weights
3. How to handle the seams between blocks?
4. The number of reads from each var seq block should be adjusted to maintain coverage


V 0.2.1 (This will leave seams)
-------------------------------

                  block 0                      block 1
    Ref seq |---------------------------|----------------------------|----------------- .....


    Var seqs
            |------------------|
            |------------------------|                 --> generate reads to give required coverage
            |------------------------------|
            |--------------|


                                        |------------------|
                                        |----------------------------------|
                                        |----------------------------|
                                        |-------------------------------------|
                                        |--------------|


Notes:
1. Depending on how annoying seams are they will be fixed in V 0.2.2


Trivia
======
Since you were dying to know: Mitty comes from James Thurber's "[The Secret Life of Walter Mitty][mitty]" one of my
favourite pieces from one of my favourite authors. Though [Wikipedia][wiki] has a less favourable interpretation of what
Walter Mitty stands for, I always thought that Thurber was celebrating the dreamer within each us who spices up the banal
parts of life with a little harmless fantasy.

[mitty]: http://www.newyorker.com/archive/1939/03/18/390318fi_fiction_thurber?currentPage=all
[wiki]: http://en.wikipedia.org/wiki/The_Secret_Life_of_Walter_Mitty