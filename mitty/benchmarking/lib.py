"""Some utility functions"""
import pandas as pd


def store_iter(fname, table_name, columns, block_size=1000000):
  """Thin wrapper around store.select

  :return: dataframe
  """
  store = pd.HDFStore(fname, 'r')
  nrows = store.get_storer(table_name).table.nrows
  for st in xrange(0, nrows, block_size):
    yield store.select(
      table_name, start=st, stop=st + block_size,
      columns=columns)
