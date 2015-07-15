#!python
from mitty.version import __version__

__cmd__ = """reads ({:s}): Generate reads from a genome (either a reference genome or a sample genome).
The read characteristics are governed by the chosen read plugin.

Commandline::

  Usage:
    reads generate <pfile>  [-v|-V] [-p]
    reads list
    reads explain (parameters|model (all|<model_name>))

  Options:
    <pfile>         Name for parameter file
    list            List available models
    explain         Explain details about the indicated model
    parameters      Explain parameters
    model           Explain a model
    all             Iterate over all models and explain them
    <model_name>    The model to explain
    -p              Show progress bar
    -v              Dump detailed logger messages
    -V              Dump very detailed logger messages
""".format(__version__)

__param__ = """Parameter file example::

  {
    "files": {
      # An absolute path is left as is
      # a relative path is taken relative to the location of the *script*
      "reference_dir": "/Users/kghose/Data/hg38/",  # Use this if the reference consists of multiple .fa files in a directory
      "reference_file": "/Users/kghose/Data/hg38/hg38.fa.gz",  # Use this if reference is a single multi-fasta file
      "dbfile": "Out/test.db"  # Genomes database file. Leave out if taking reads from a reference, or from VCF file
      "input_vcf": "Out/test.vcf",  # Use this if using a VCF file. Leave out if taking reads from a reference, or from VCF file
                                    # The vcf file is converted to a Mitty genome database that is stored alognside the
                                    # fastq files
      "output_prefix": "Out/reads", # Output file name prefix
                                    # the reads will be called reads.fq and reads_c.fq if we call for corrupted reads too
      "interleaved": true     # Set to false if you need separate files for each mate of a pair
                              # files will then be named reads_1.fq, reads_c_1.fq and reads_2.fq, reads_c_2.fq
                              # If the reads are actually not paired and you set interleaved to false you will get two
                              # files _1 and _2 and all the data will be in _1 only
    },
    "sample_name": "g0_s0",   # Name of sample
    "rng": {
      "master_seed": 1
    },
    "chromosomes": [1,2],               # List the chromosomes the reads should be taken from
    "variants_only": false,             # If true, reads will only come from the vicinity of variants
    "variant_window": 500,              # Only used if variants_only is true. Reads will be taken from this vicinity
    "corrupt": true,                    # If true, corrupted reads will also be written
    "coverage": 10,                     # Coverage
    "coverage_per_block": 0.01          # For each block of the simulation we generate reads giving this coverage
                                        # In some simulations we will have too many reads to fit in memory and we must
                                        # write them out in blocks of smaller size.
    "read_model": "simple_illumina",    # Model specific parameters, need to be under the key "model_params"
    "model_params": {
      "read_len": 100,          # length of each read
      "template_len_mean": 250, # mean length of template
      "template_len_sd": 30,    # sd of template length
      "max_p_error": 0.01,      # Maximum error rate at tip of read
      "k": 20                   # exponent
    }

  }
"""
__doc__ = __cmd__ + __param__
# We split this up so that when we print the help, we can print it as separate pages. It got cluttered when printed all
# together. We want this to be the docstring because that is a nice format for our generated documentation

from contextlib import contextmanager
import json
import os
import time

import numpy as np
import docopt
import click

import mitty.lib
import mitty.lib.reads as lib_reads
import mitty.lib.io as mio
import mitty.lib.variants as vr

import logging
logger = logging.getLogger(__name__)


class ReadSimulator:
  """A convenience class that wraps the parameters and settings for a read simulation"""
  def __init__(self, base_dir, params):
    """Create a read simulator object

    :param base_dir: the directory with respect to which relative file paths will be resolved
    :param params: dict loaded from json file
    """

    fname_prefix = mitty.lib.rpath(base_dir, params['files']['output_prefix'])
    if not os.path.exists(os.path.dirname(fname_prefix)):
      os.makedirs(os.path.dirname(fname_prefix))

    self.ref = mio.Fasta(multi_fasta=mitty.lib.rpath(base_dir, params['files'].get('reference_file', None)),
                         multi_dir=mitty.lib.rpath(base_dir, params['files'].get('reference_dir', None)),
                         persistent=True)

    self.sample_name = params.get('sample_name', None)
    if 'dbfile' in params['files']:
      pop_db_name = mitty.lib.rpath(base_dir, params['files']['dbfile'])
      self.pop = vr.Population(fname=pop_db_name)
    else:
      self.pop = None
      logger.debug('Taking reads from reference')

    master_seed = int(params['rng']['master_seed'])
    assert 0 < master_seed < mitty.lib.SEED_MAX
    self.seed_rng = np.random.RandomState(seed=master_seed)

    self.chromosomes = params['chromosomes']
    self.coverage = float(params['coverage'])

    total_blocks_to_do = self.coverage / float(params['coverage_per_block'])

    chrom_meta = self.ref.get_seq_metadata()
    sum_of_chromosome_lengths = float(sum([chrom_meta[chrom - 1]['seq_len'] for chrom in self.chromosomes]))
    self.blocks_for_chromosome = {chrom: int(max(1, round(total_blocks_to_do * chrom_meta[chrom - 1]['seq_len'] / sum_of_chromosome_lengths)))
                                  for chrom in self.chromosomes}

    self.read_model = mitty.lib.load_reads_plugin(params['read_model']).Model(**params['model_params'])
    self.corrupt_reads = bool(params['corrupt'])

    self.variants_only = params.get('variants_only', None)
    self.variant_window = int(params.get('variant_window', 200)) if self.variants_only else None

    if params['files'].get('interleaved', True):
      self.fastq_fp = [open(fname_prefix + '.fq', 'w')] * 2
      self.fastq_c_fp = [open(fname_prefix + '_c.fq', 'w')] * 2 if self.corrupt_reads else [None, None]
    else:  # Need two files separately
      self.fastq_fp = [open(fname_prefix + '_1.fq', 'w'), open(fname_prefix + '_2.fq', 'w')]
      self.fastq_c_fp = [open(fname_prefix + '_c_1.fq', 'w'), open(fname_prefix + '_c_2.fq', 'w')] if self.corrupt_reads else [None, None]

    self.first_serial_no = 0

  def get_total_blocks_to_do(self):
    return sum(self.blocks_for_chromosome.values()) * 2  # Two copies for each chromosome

  def get_chromosome_list(self):
    return self.chromosomes

  def get_blocks_to_do(self, chrom):
    return self.blocks_for_chromosome

  def generate_and_save_reads(self, chrom, cpy):
    """Grab the appropriate seq, apply variants as needed, generate reads, roll cigars and then save."""
    if self.pop is not None:
      ml = self.pop.get_master_list(chrom=chrom)
      v_index = self.pop.get_sample_chromosome(chrom=chrom, sample_name=self.sample_name)
    else:
      ml, v_index = vr.VariantList(), []  # Need a dummy variant list for nulls

    seq, variant_waypoints = lib_reads.expand_sequence(self.ref[chrom]['seq'], ml, v_index, cpy)
    seq_c = mitty.lib.string.translate(seq, mitty.lib.DNA_complement)
    coverage_per_block = 0.5 * self.coverage / self.blocks_for_chromosome[chrom]
    for blk in range(self.blocks_for_chromosome[chrom]):
      reads, paired = generate_reads(seq, seq_c, variant_waypoints,
                                     self.variant_window, self.read_model,
                                     coverage_per_block,
                                     self.corrupt_reads,
                                     self.seed_rng,
                                     variants_only=self.variants_only)
      pos, cigars = lib_reads.roll_cigars(variant_waypoints, reads)
      template_count = write_reads_to_file(self.fastq_fp[0], self.fastq_fp[1],
                                           self.fastq_c_fp[0], self.fastq_c_fp[1],
                                           reads, paired, pos, cigars,
                                           chrom, cpy,
                                           self.first_serial_no)
      self.first_serial_no += template_count
      yield


def generate_reads(seq, seq_c, variant_waypoints, variant_window,
                   read_model, coverage, corrupt, seed_rng, variants_only=False):
  """Wrapper around read function to handle both regular reads as well as reads restricted to around

  :param seq:      forward sequence
  :param seq_c:    complement sequence
  :param variant_waypoints: as returned by expand_sequence
  :param variant_window: how many bases before and after variant should we include
  :param corrupt:  T/F generate corrupted reads too or not
  :param seed_rng:    rng for seed generation
  :param variants_only: set True if we want reads only from variant regions
  :return:
  """
  if variants_only:
    reads, paired = [], False
    for v in variant_waypoints:  # [1:-1]:
      start, stop = max(v[1] - variant_window, 0), min(v[1] + variant_window, len(seq))
      # v[1] is the pos of the variant in sequence coordinates (rather than ref coordinates)
      these_reads, paired = read_model.get_reads(seq, seq_c,
                                                 start_base=start, end_base=stop,
                                                 coverage=coverage,
                                                 corrupt=corrupt,
                                                 seed=seed_rng.randint(0, mitty.lib.SEED_MAX))
      reads += [these_reads]
    reads = np.concatenate(reads)
  else:
    reads, paired = read_model.get_reads(seq, seq_c,
                                         coverage=coverage,
                                         corrupt=corrupt,
                                         seed=seed_rng.randint(0, mitty.lib.SEED_MAX))
  return reads, paired


def write_reads_to_file(fastq_fp_1, fastq_fp_2,
                        fastq_c_fp_1, fastq_c_fp_2,
                        reads, paired, pos, cigars,
                        chrom, cpy,
                        first_serial_no):
  """
  :param fastq_fp_1:   file pointer to perfect reads file
  :param fastq_fp_2:   file pointer to perfect reads file 2. Same as 1 if interleaving. None if not paired
  :param fastq_c_fp_1: file pointer to corrupted reads file. None if no corrupted reads being written
  :param fastq_c_fp_2: file pointer to corrupted reads file 2. Same as 1 if interleaving. None if no corrupted reads being written. None if not paired
  :param reads:        reads recarray
  :param paired:       bool, are reads paired
  :param pos:          list of POS values for the reads
  :param cigars:       list of cigar values for the reads
  :param chrom:           chromosome number
  :param cpy:           chromosome copy
  :param first_serial_no: serial number of first template in this batch
  :return: next_serial_no: the serial number the next batch should start at
  """
  #cntr = first_serial_no
  ro, pr_seq, cr_seq, phred = reads['read_order'], reads['perfect_reads'], reads['corrupt_reads'], reads['phred']

  if paired:
    for n in xrange(0, reads.shape[0], 2):
      qname = '{:d}|{:d}|{:d}|{:d}|{:d}|{:s}|{:d}|{:d}|{:s}'.format(first_serial_no + n/2, chrom, cpy, ro[n], pos[n], cigars[n], ro[n + 1], pos[n + 1], cigars[n + 1])
      fastq_fp_1.write('@' + qname + '/1\n' + pr_seq[n] + '\n+\n' + '~' * len(pr_seq[n]) + '\n')
      fastq_fp_2.write('@' + qname + '/2\n' + pr_seq[n + 1] + '\n+\n' + '~' * len(pr_seq[n + 1]) + '\n')
      if fastq_c_fp_1 is not None:
        fastq_c_fp_1.write('@' + qname + '/1\n' + cr_seq[n] + '\n+\n' + phred[n] + '\n')
        fastq_c_fp_2.write('@' + qname + '/2\n' + cr_seq[n + 1] + '\n+\n' + phred[n + 1] + '\n')
    template_count = reads.shape[0] / 2
  else:
    for n in xrange(0, reads.shape[0]):
      qname = '{:d}|{:d}|{:d}|{:d}|{:d}|{:s}'.format(first_serial_no + n, chrom, cpy, ro[n], pos[n], cigars[n])
      fastq_fp_1.write('@' + qname + '\n' + pr_seq[n] + '\n+\n' + '~' * len(pr_seq[n]) + '\n')
      if fastq_c_fp_1 is not None:
        fastq_c_fp_1.write('@' + qname + '\n' + cr_seq[n] + '\n+\n' + phred[n] + '\n')
    template_count = reads.shape[0]
  return template_count


@contextmanager
def no_bar(**kwargs):
  class NoBar():
    def update(self, _):
      pass
  yield NoBar()


@click.command()
@click.argument('param_fname', type=click.Path(exists=True))
@click.option('-v', count=True, help='Verbosity level')
@click.option('-p', is_flag=True, help='Show progress bar')
def cli(param_fname, v, p):
  """Generate reads (fastq) given a parameter file"""
  level = logging.DEBUG if v > 1 else logging.WARNING
  logging.basicConfig(level=level)
  if v == 1:
    logger.setLevel(logging.DEBUG)

  base_dir = os.path.dirname(param_fname)     # Other files will be with respect to this
  params = json.load(open(param_fname, 'r'))

  simulation = ReadSimulator(base_dir, params)

  bar_func = click.progressbar if p else no_bar
  t0 = time.time()
  with bar_func(length=simulation.get_total_blocks_to_do(), label='Generating reads') as bar:
    for chrom in simulation.get_chromosome_list():
      for cpy in [0, 1]:
        for _ in simulation.generate_and_save_reads(chrom, cpy):
          bar.update(1)
  t1 = time.time()
  logger.debug('Took {:f}s to write {:d} templates'.format(t1 - t0, simulation.first_serial_no))


def print_list(cmd_args):
  print('\nAvailable models\n----------------')
  for name, mod_name in mitty.lib.discover_all_reads_plugins():
    print('- {:s} ({:s})\n'.format(name, mod_name))


def explain(cmd_args):
  if cmd_args['parameters']:
    print(__param__)
  elif cmd_args['model']:
    explain_read_model(cmd_args['<model_name>'])


def explain_read_model(name):
  try:
    mod = mitty.lib.load_reads_plugin(name)
  except ImportError:
    print('No model named {:s}'.format(name))
    return
  try:
    print(mod._description)
    print(mitty.lib.model_init_signature_string(mod.Model.__init__))
  except AttributeError:
    print('No help for model "{:s}" available'.format(name))
  return


def explain_all_read_models():
  for name, mod_name in mitty.lib.discover_all_reads_plugins():
    explain_read_model(name)


if __name__ == '__main__':
  cli()