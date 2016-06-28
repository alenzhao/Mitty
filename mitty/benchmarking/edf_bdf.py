"""Load in bdf and evdf and find reads under the calls


python ~/Code/Mitty/mitty/benchmarking/evdf-bdf.py eval.df.csv.gz test-bam.csv.gz calls-reads.csv -v


"""
import logging
import time
from multiprocessing import Pool
import os
import gzip

import click
import pandas as pd

from mitty.benchmarking.dfcols import *

import numpy
import pyximport
pyximport.install(setup_args={"include_dirs": numpy.get_include()})
from mitty.benchmarking.edf_bdf_cy import get_templates_over_calls, get_templates_over_calls_old


logger = logging.getLogger(__name__)


@click.command()
@click.argument('evdf')
@click.argument('bdf')
@click.argument('outcsv')
@click.option('-t', default=4, help='Threads')
@click.option('--block-size', default=100000, help='Block size')
@click.option('-v', count=True, help='Verbosity level')
def cli(evdf, bdf, outcsv, t, block_size, v):
  """Associate reads in a BAM with variant calls."""
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  if os.path.exists(outcsv):
    # logger.error('Output file {} already exists. Will not over-write'.format(outcsv))
    # exit(1)
    logger.warning('Output file {} already exists, removing'.format(outcsv))
    os.remove(outcsv)

  unsorted_file = 'unsorted-' + outcsv
  process(evdf, bdf, unsorted_file, t, block_size)
  sort_df(unsorted_file, outcsv)


def process(evdf, bdf, outcsv, t, block_size):
  eval_df = pd.read_csv(evdf, dtype={'chrom': 'S10'}, compression='gzip' if evdf.endswith('gz') else None)
  logger.debug('Loaded {}'.format(evdf))
  g = ((bam_df, eval_df) for bam_df in pd.read_csv(bdf, compression='gzip' if bdf.endswith('gz') else None, chunksize=block_size))

  t0 = time.time()
  p = Pool(t)
  write_header(outcsv)
  total_reads, total_lines = 0, 0
  for l in p.imap_unordered(get_all_templates_over_calls, g):
    pd.DataFrame(l, columns=calls_reads_cols).to_csv(
      outcsv, index=False, header=False, compression='gzip' if outcsv.endswith('gz') else None, mode='a')
    total_reads += block_size
    total_lines += len(l)
    logger.debug('{} lines'.format(total_lines))
  t1 = time.time()
  logger.debug('Took {:0.10}s to process {}/{} templates/calls ({:0.10} templates/s) with {} threads'.
               format(t1 - t0, total_reads, len(eval_df), total_reads/(t1 - t0), t))


def sort_df(incsv, outcsv, remove=True):
  """

  :param incsv:  Unsorted data frame
  :param outcsv: Sorted data frame
  :param remove: Remove unsorted file when done
  :return:
  """
  logger.debug('Loading unsorted file and sorting')
  pd.read_csv(incsv).sort_values('call_p1').to_csv(
    outcsv, index=False, compression='gzip' if outcsv.endswith('gz') else None)
  logger.debug('Saved sorted data frame to {}'.format(outcsv))
  if remove:
    os.remove(incsv)
    logger.debug('Removing unsorted file {}'.format(incsv))


def write_header(outcsv):
  if outcsv.endswith('gz'):
    with gzip.open(outcsv, 'w') as fp:
      fp.write(','.join(calls_reads_cols) + '\n')
  else:
    with open(outcsv, 'w') as fp:
      fp.write(','.join(calls_reads_cols) + '\n')


def get_all_templates_over_calls(args):
  bam_df, eval_df = args
  import traceback
  try:
    return get_templates_over_calls(bam_df, eval_df)
  except Exception as e:
    traceback.print_exc()
    print()
    raise e


if __name__ == '__main__':
  cli()