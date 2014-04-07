"""Example parameter file for reads program

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""

model_name = 'perfect'  # Refers to the stock reads model Plugins/Reads/perfect_plugin.py
#seq_header = '11 dna:chromosome chromosome:GRCh37:11:1:135006516:1'
#seq_file = 'Data/mutated_human_chrom_11.smalla'
seq_header = 'gi|4630864|dbj|AB026117.1| Porcine adenovirus 3 DNA, complete genome'
seq_file = 'Data/porcine_circovirus.smalla'
corrupted_reads_file = 'Data/corrupted_reads.bam'
perfect_reads_file = 'Data/perfect_reads.bam'
output_type = 'bam'
start = 0  # We can restrict the region the reads are taken from
stop = -1
coverage = 5
average_read_len = 100  # This is used by reads.py to figure out how many reads to generate for given coverage
args = {
  'paired': False, #True,
  'read_len': 100,
  'template_len': 1000,
  'rng_seed': 1            # Only one RNG
}