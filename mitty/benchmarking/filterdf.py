"""MQ and more"""
import logging
import time

import click
import numpy as np
import pandas as pd

import mitty.benchmarking.lib as mlib


logger = logging.getLogger(__name__)


@click.command()
@click.version_option()
@click.argument('df_fname', type=click.Path(exists=True))
@click.argument('filter', type=str)
@click.argument('outdf', type=click.Path())
@click.option('--block-size', default=1000000, help='Block size')
@click.option('-v', count=True, help='Verbosity level')
def cli(df_fname, filter, outdf, block_size, v):
  """Given a filter string, apply it to the df and push out a copy of the df with that filter applied"""
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  st = pd.HDFStore(df_fname, 'r')
  tn = 'bdf' if 'bdf' in st else 'crdf'
  data_columns = [c for c, v in st.get_storer(tn).table.colindexed.items() if v and c != 'index']
  nrows = st.get_storer(tn).table.nrows
  st.close()

  st = pd.HDFStore(outdf, mode='a', complevel=9, complib="blosc", format='t')
  t0 = time.time()
  for df in mlib.store_iter(
    df_fname, tn, columns=None, block_size=block_size):
    st.append(tn, df.query(filter), data_columns=data_columns, index=False)
    t1 = time.time()
    logger.debug('Took {:0.10}s to process {} rows ({:0.10} rows/s)'.format(t1 - t0, df.shape[0], df.shape[0] / (t1 - t0)))

  t1 = time.time()
  logger.debug('Took {:0.10}s to process {} rows ({:0.10} rows/s)'.format(t1 - t0, nrows, nrows / (t1 - t0)))

  st.create_table_index(tn, optlevel=9, kind='full')
  t2 = time.time()
  logger.debug('Took {:0.10}s to index'.format(t2 - t1))


if __name__ == '__main__':
  cli()