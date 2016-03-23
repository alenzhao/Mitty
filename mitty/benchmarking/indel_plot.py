from itertools import izip
import time
import json
import io

import click
import numpy as np
import pysam

import mitty.lib.variants as vr
from mitty.benchmarking.creed import MAX_D_ERROR

import logging
logger = logging.getLogger(__name__)


@click.group()
@click.version_option()
def cli():
  """Plot alignment accuracy categorized by variant type."""
  pass


@cli.command()
@click.argument('perbam', type=click.Path(exists=True))
@click.argument('indeljson', type=click.Path())
@click.option('--indel-range', help='Maximum base pair count of indels we process', type=int, default=50)
@click.option('-p', is_flag=True, help='Show progress bar')
@click.option('-v', count=True, help='Verbosity level')
def plot(perbam, indeljson, indel_range, p, v):
  """INDELBAM -> plot + json:

  \b
  Given a perfect BAM file with indel annotations added by indelbam, use the information in the
  extended tags to work out alignment accuracy parametrized by the sample indel lengths. Generate
  a plot and a json file."""
  perbam_fp = pysam.AlignmentFile(perbam, 'rb')
  total_read_count = perbam_fp.mapped + perbam_fp.unmapped  # Sadly, this is only approximate
  progress_bar_update_interval = int(0.01 * total_read_count)

  t0 = time.time()
  with click.progressbar(
    length=total_read_count, label='Counting reads', file=None if p else io.BytesIO()) as bar:
    json.dump(
      parse_bam(perbam_fp, indel_range, bar.update, progress_bar_update_interval),
      open(indeljson, 'w'), cls=NumpyJsonEncoder)
  t1 = time.time()
  logger.debug('Processed {:d} reads in {:2.2f}s'.format(total_read_count, t1 - t0))


@cli.command()
@click.argument('indeljson', type=click.Path(exists=True))
@click.option('--indel-range', help='Maximum base pair count of indels we process', type=int, default=50)
def replot(indeljson):
  """json -> plot:

  \b
  Given a json file with alignment accuracy parametrized by the sample indel lengths,
  generate a plot."""
  print(json)


def parse_bam(perbam_fp, indel_range, progress_callback=None, progress_bar_update_interval=None):
  """Given a BAM, analyze the extended tags and score the read accuracy, parametrized by indel length."""
  variants = np.zeros((2 * indel_range + 1,), dtype=[('variant_length', 'i2'),
                                                     ('novel_variant_count', 'i4'),
                                                     ('known_variant_count', 'i4')])
  variants['variant_length'] = np.arange(-indel_range, indel_range + 1, dtype=int)

  # Dimensions are: indel_size, d_error, read_count
  reference_read_counts = np.zeros((MAX_D_ERROR + 1), dtype=int)
  novel_variant_read_counts = np.zeros((2 * indel_range + 1, MAX_D_ERROR + 1), dtype=int)
  known_variant_read_counts = np.zeros((2 * indel_range + 1, MAX_D_ERROR + 1), dtype=int)

  cnt0, cnt = 0, 0
  for cnt, r in enumerate(perbam_fp):
    d_e = r.get_tag('Xd')
    #corr = r.get_tag('Xf')
    if r.has_tag('Z0'):  # read covers at least one variant
      for vs, vf in zip(r.get_tag('Z0'), r.get_tag('Z1')):
        if -indel_range <= vs <= indel_range:
          idx = vs + indel_range
          if vf == 1:
            variants['novel_variant_count'][idx] += 1
            novel_variant_read_counts[idx, d_e] += 1
          elif vf == 2:
            variants['known_variant_count'][idx] += 1
            known_variant_read_counts[idx, d_e] += 1

    else:
      reference_read_counts[d_e] += 1

    if (cnt + 1) % progress_bar_update_interval:
      progress_callback(cnt - cnt0)
      cnt0 = cnt

  return {
    'variants': variants,
    'reference_read_counts': reference_read_counts,
    'novel_variant_read_counts': novel_variant_read_counts,
    'known_variant_read_counts': known_variant_read_counts
  }


class NumpyJsonEncoder(json.JSONEncoder):
  def default(self, obj):
    if isinstance(obj, np.ndarray):
      return obj.tolist()
    return json.JSONEncoder.default(self, obj)






# @cli.command()
# @click.option('--perbam', type=click.Path(exists=True), help='PerBAM file annotated by indelbam')
# @click.option('--indel-json', type=click.Path(), help='indel-json file')
#
# @click.argument('perbam', type=click.Path(exists=True))
# @click.argument('vdb', type=click.Path(exists=True))  # , help='File name of genome database')
# @click.argument('out-json', type=click.Path())
# @click.option('--sample-name', help='Name of sample to compare against. Leave out to test against population')
# @click.option('--indel-range', help='Maximum base pair count of indels we process', type=int, default=50)
# @click.option('-p', is_flag=True, help='Show progress bar')
# def cli(perbam, out_json, vdb, sample_name, indel_range, p):
#
#   bam_fp = pysam.AlignmentFile(perbam, 'rb')
#   pop = vr.Population(vdb)
#   cat_read_counts = None
#   with click.progressbar(length=2 * len(pop.get_chromosome_list()), label='Analyzing indels', file=None if p else io.BytesIO()) as bar:
#     for ch in pop.get_chromosome_list():
#       # cat_read_counts = categorize_data_from_one_chromosome(
#       #   bam_fp, pop, ch, sample_name=sample_name, cat_read_counts=cat_read_counts, max_indel=indel_range)
#       for cat_read_counts in categorize_data_from_one_chromosome(
#           bam_fp, pop, ch, sample_name=sample_name, cat_read_counts=cat_read_counts, max_indel=indel_range):
#         bar.update(1)
#
#   json.dump(cat_read_counts, open(out_json, 'w'), cls=NumpyJsonEncoder)


if __name__ == '__main__':
  cli()