Ramblings AKA. Diary of a Developer
===================================

(In which I write down some of my design decisions.)

Simulating a million Genomes in Python
--------------------------------------

Memory considerations
+++++++++++++++++++++

In the early days Python was a no-brainer in terms of prototyping ideas. It got us off the ground very fast, it allowed
me to write comprehensive and compact code and to use things like docopts and sphinx that allowed good command line
interfaces and painless, excellent documentation generation.

The real challenge started when we were thinking of population simulations. The primary consideration was that, being a
simulation meant for benchmarking and tool testing, the code should be accurate, rather than fast, and if we wanted to
simulate a million genomes, we should be able to do it i.e. we would sacrifice time for space.

First, some numbers:

A single variant requires at least 7 bytes.

pos      4 bytes
ref      1 byte
alt      1 byte
genotype 1 byte

For speed of computation we like to keep a stop field which says where the variant stops (mostly pos + len(ref)) bringing
us to 11 bytes

pos      4 bytes
stop     4 bytes
ref      1 byte
alt      1 byte
genotype 1 byte

A single genome has about 3e6 (three million) variants (you see we are only considering snps, but this gets us a good
estimate of magnitude). This means that a single genome requires between 21e6 and 33e6 bytes.

The 1000G data set has about 100e6 (hundred million) unique variants from 2500 samples, and a million genome data will
probably have a billion.

Considering a population of only 2500 we see that

If we store each genome separately, we will require between 52.5 and 82.5 billion bytes. If we store just the unique
variants we will end up requiring 0.7 to 1.1 billion bytes + something for the genotypes of each sample.

How would the second scheme work?

We would use 6 or 10 bytes to store just the variant (no genotype information) and then for each sample just store the
genotype information and a pointer to the variant, which takes the form:

index    4 bytes  (points to the variant in the master list)
genotype 1 byte

(In practice, I compressed this down to 4 bytes, by using 30 bits for the index and 2 bits for the diploid genotype)

0.6 to 1.0 billion bytes for the master list of variants.
15 million bytes per sample = 37.5 billion bytes for 2500 samples.

So it is clear that even with a "tiny" sample like the 1000G data we can not fit the whole population in memory at the
same time, but by using a master list of unique variants we can halve our memory requirements.


Sampling a population into existence
++++++++++++++++++++++++++++++++++++

While I was


Report generation
-----------------
Static textual html pages are pretty easy to generate with pure Python. HTML pages with static or interactive graphs are a
little harder but some distance can be covered using Matplotlib_ and the mpld3_ library. If we are willing to abandon
Matplotlib, we can use Bokeh_ which is perfect for this. Matplotlib, however, has a high bar for publication type static
figures that are written to pdf.

The challenge comes when we
need interactive plots that show large datasets. Loading all the data in at once causes the page to become sluiggish.

In a pure Matplotlib interactive figure we can use callbacks to load views of the data from disk as needed. There is no
way to translate this to a web-page, currently.

Using Bokeh we can build such a plot, but we need a bokeh server to serve the page and
respond to the interaction commands. This is not feasible with our current platform.

From these considerations I decided to use the following approach for the platform:

1. There is a single page html report generated that can be easily seen on the platform
2. The html has some interactive textual tables summarizing benchmark data.
3. There are a few interactive summary plots which attempt to give a graphical overview of the errors. The summary plots
bin the data so that we handle a small quantity of numbers which keep our web-page snappy.

For people using the terminal, we generate interactive plots in Bokeh/Matplotlib that load data from disk as needed to
display detailed views of the errors for extensive debugging.

.. _Matplotlib: http://matplotlib.com
.. _mpld3: http://mpld3.com
.. _Bokeh: http://bokeh.io
