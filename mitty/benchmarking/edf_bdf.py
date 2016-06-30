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
from mitty.benchmarking.edf_bdf_cy import get_templates_over_calls


logger = logging.getLogger(__name__)


@click.command()
@click.argument('evdf')
@click.argument('bdf')
@click.argument('outh5')
@click.option('-t', default=4, help='Threads')
@click.option('--block-size', default=100000, help='Block size')
@click.option('--max-blocks-to-do', type=int, help='Maximum blocks to do (for debugging, single thread only)')
@click.option('-v', count=True, help='Verbosity level')
def cli(evdf, bdf, outh5, t, block_size, max_blocks_to_do, v):
  """Associate reads in a BAM with variant calls."""
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  if os.path.exists(outh5):
    # logger.error('Output file {} already exists. Will not over-write'.format(outcsv))
    # exit(1)
    logger.warning('Output file {} already exists, removing'.format(outh5))
    os.remove(outh5)

  process(evdf, bdf, outh5, t, block_size, max_blocks_to_do)
  # sort_df(unsorted_file, outcsv)


# TODO: perhaps have a lockable HDF5 store in the threads?
def process(evdf, bdf, outh5, t, block_size, max_blocks_to_do=None):
  # TODO, turn this into HDF5 too
  eval_df = pd.read_csv(evdf, dtype={'call_chrom': 'S10'}, compression='gzip' if evdf.endswith('gz') else None)
  logger.debug('Loaded {}'.format(evdf))
  g = ((bam_df, eval_df) for bam_df in pd.read_csv(bdf, compression='gzip' if bdf.endswith('gz') else None, chunksize=block_size))

  t0 = time.time()
  p = Pool(t)
  st = pd.HDFStore(outh5, mode='a', complevel=9, complib="blosc", format='t')

  total_reads, total_lines, blocks = 0, 0, 0
  for l in p.imap_unordered(get_all_templates_over_calls, g):
    # NOTES:
    # strings and minitemsize: http://stackoverflow.com/questions/15988871/hdfstore-appendstring-dataframe-fails-when-string-column-contents-are-longer
    # http://pandas-docs.github.io/pandas-docs-travis/io.html#string-columns

    # This commented code is the .to_hdf version. It works fine, but may be creating the index each time
    # pd.DataFrame(l, columns=calls_reads_cols).to_hdf(
    #   outh5, 'edfbdf', format='t', append=True,
    #   min_itemsize=1000
    # )
    st.append('edfbdf',
              pd.DataFrame(l, columns=edf_bdf_call_cols + edf_bdf_read_cols),
              data_columns=edf_bdf_data_cols, index=False)

    total_reads += block_size
    total_lines += len(l)
    logger.debug('{} lines'.format(total_lines))

    blocks += 1
    if max_blocks_to_do is not None and blocks >= max_blocks_to_do:
      break

  logger.debug('Creating index')
  # http://pandas-docs.github.io/pandas-docs-travis/io.html#indexing
  st.create_table_index('edfbdf', optlevel=9, kind='full', columns=edf_bdf_data_cols)

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