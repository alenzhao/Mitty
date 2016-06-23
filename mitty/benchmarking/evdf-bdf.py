"""Load in bdf and evdf and find reads under the calls


python ~/Code/Mitty/mitty/benchmarking/evdf-bdf.py eval.df.csv.gz test-bam.csv.gz calls-reads.csv -v


"""
import logging
import time
from multiprocessing import Pool
import os

import click
import pandas as pd

from mitty.benchmarking.dfcols import *

import numpy
import pyximport
pyximport.install(setup_args={"include_dirs": numpy.get_include()})
from mitty.benchmarking.eval_bdf_cy import *


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
    logger.error('Output file {} already exists. Will not over-write'.format(outcsv))
    exit(1)

  eval_df = pd.read_csv(evdf, dtype={'chrom': 'S10'}, compression='gzip' if evdf.endswith('gz') else None)
  g = ((bam_df, eval_df) for bam_df in pd.read_csv(bdf, compression='gzip' if bdf.endswith('gz') else None, chunksize=block_size))

  t0 = time.time()
  p = Pool(t)
  write_header(outcsv)
  total_reads = 0
  for l in p.imap_unordered(get_all_templates_over_calls, g):
    pd.DataFrame(l, columns=calls_reads_cols).to_csv(
      outcsv, index=False, header=False, compression='gzip' if outcsv.endswith('gz') else None, mode='a')
    total_reads += block_size
  t1 = time.time()
  logger.debug('Took {:0.10}s to process {}/{} templates/calls ({:0.10} templates/s) with {} threads'.
               format(t1 - t0, total_reads, len(eval_df), total_reads/(t1 - t0), t))


def write_header(outcsv):
  with open(outcsv, 'w') as fp:
    fp.write(','.join(calls_reads_cols) + '\n')


def get_all_templates_over_calls(args):
  bam_df, eval_df = args
  call_read_rows = []
  for mate in ['m1', 'm2']:
    for mode in ['a', 'c']:
      call_read_rows += get_templates_over_calls(bam_df, eval_df, mate, mode)
  return call_read_rows


# DEBUG:__main__:Took 19.30675101s to process 60000/244546 reads/calls (3107.72123 reads/s) with 2 threads
def get_templates_over_calls_python(bam_df, eval_df, mate='m1', mode='a'):
  """

  :param bam_df: Whole or part of a bdf
  :param eval_df: Whole (or part) of an evdf
  :param mate: m1 or m2 : sort on mate1 or mate2
  :param mode: 'a' or 'c' : sort on aligned or correct position
  :return:
  """
  t0 = time.time()
  call_read_rows = []

  feature_key = [('{}_{}_under_feature'.format(mate, mode), 1)]

  p1_key = '{}_{}_p1'.format(mate, mode)  # e.g. 'm1_a_p1'
  p2_key = '{}_{}_p2'.format(mate, mode)  # e.g. 'm1_a_p2'

  d_sorted = bam_df.sort_values(by=[p1_key])

  r_p1 = d_sorted[p1_key].values
  r_p2 = d_sorted[p2_key].values

  v_p1 = eval_df['call_p1'].values
  v_p2 = eval_df['call_p2'].values

  call_cnt = len(eval_df)
  t_idx1, t_idx2, template_cnt = 0, 0, len(bam_df)

  for call_index in range(call_cnt):
    if t_idx1 >= template_cnt:
      break

    # Move read window such that all reads in window are over call
    while t_idx1 < template_cnt and (r_p2[t_idx1] < v_p1[call_index]):  # Move upper edge of window
      t_idx1 += 1

    while t_idx2 < template_cnt and (r_p1[t_idx2] < v_p2[call_index]):  # Move lower edge of window
      t_idx2 += 1
      # There are some shenanigans here with the last read in the file
      # (we don't consider it) but we don't care.

    if t_idx2 > t_idx1:
      citems = eval_df.iloc[call_index].to_dict().items()
      for r in d_sorted.iloc[t_idx1:t_idx2].iterrows():
        call_read_rows.append(
          {k: v for k, v in citems + r[1].to_dict().items() + feature_key}
        )

  t1 = time.time()
  logger.debug('Matched {} reads against {} calls in {:0.3} s'.format(template_cnt, call_cnt, t1 - t0))

  return call_read_rows


if __name__ == '__main__':
  cli()