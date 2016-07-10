"""Load in bdf and evdf and find reads under the calls
python ~/Code/Mitty/mitty/benchmarking/evdf-bdf.py eval.df.csv.gz test-bam.csv.gz calls-reads.csv -v
"""
import logging
import time
from multiprocessing import Pool
import os

import click
import pandas as pd

import mitty.benchmarking.dfcols as dfcols

import numpy
import pyximport
pyximport.install(setup_args={"include_dirs": numpy.get_include()})
from mitty.benchmarking.edf_bdf_cy import get_templates_over_calls


logger = logging.getLogger(__name__)


@click.command()
@click.version_option()
@click.argument('evdf')
@click.argument('bdf')
@click.argument('outh5')
@click.option('-t', default=4, help='Threads')
@click.option('--block-size', default=1000000, help='Block size')
@click.option('--max-templates', type=int, help='Maximum templates to do (for debugging, quantized to "block_size" increments)')
@click.option('--call-hash', type=int, help='Compute edf-bdf for this call only (for debugging)')
@click.option('-v', count=True, help='Verbosity level')
def cli(evdf, bdf, outh5, t, block_size, max_templates, call_hash, v):
  """Associate reads in a BAM with variant calls."""
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  if os.path.exists(outh5):
    # logger.error('Output file {} already exists. Will not over-write'.format(outcsv))
    # exit(1)
    logger.warning('Output file {} already exists, removing'.format(outh5))
    os.remove(outh5)

  process(evdf, bdf, outh5, t, block_size, max_templates, call_hash)


def bdf_iter(edf, bdf_st, block_size, max_templates=None):
  nrows = bdf_st.get_storer('bdf').table.nrows

  for offset in xrange(
    0,
    nrows if max_templates is None else min(nrows, max_templates),
    block_size):
    yield {
      'bdf': bdf_st.select('bdf', start=offset, stop=offset + block_size),
      'edf': edf
    }


def process(edf_name, bdf_name, outh5, t, block_size, max_templates=None, call_hash=None):

  edf = pd.read_hdf(edf_name, 'edf', columns=dfcols.get_edf_data_cols().keys())  # We only need these ones here
  logger.debug('Loaded {}'.format(edf_name))
  if call_hash is not None:
    edf = edf[edf['call_hash'] == call_hash]
    logger.warning('Only processing call with call_hash {}'.format(call_hash))

  bdf_st = pd.HDFStore(bdf_name, mode='r')
  nrows = bdf_st.get_storer('bdf').table.nrows

  t0 = time.time()
  p = Pool(t)
  st_out = pd.HDFStore(outh5, mode='a', complevel=9, complib="blosc", format='t')

  total_lines = 0
  for cr in p.imap_unordered(get_all_templates_over_calls, bdf_iter(edf, bdf_st, block_size, max_templates)):
    if cr is not None:
      st_out.append('crdf', pd.DataFrame(cr),
                    data_columns=dfcols.get_edf_data_cols().keys(), index=False)

      total_lines += len(cr)
    logger.debug('{} lines'.format(total_lines))

  logger.debug('Creating index')
  # http://pandas-docs.github.io/pandas-docs-travis/io.html#indexing
  st_out.create_table_index('crdf', optlevel=9, kind='full')

  t1 = time.time()
  logger.debug('Took {:0.10}s to process {}/{} templates/calls ({:0.10} templates/s) with {} threads'.
               format(t1 - t0, nrows, len(edf), nrows / (t1 - t0), t))

  # Store header (with sequence meta info) here, so we can write out a VCF with the call information
  # as needed for intersection/difference.


def get_all_templates_over_calls(args):
  import traceback
  try:
    return get_templates_over_calls(args['bdf'], args['edf'])
  except Exception as e:
    traceback.print_exc()
    print()
    raise e


if __name__ == '__main__':
  cli()