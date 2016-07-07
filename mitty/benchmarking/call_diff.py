import click
import logging
import time
import json
import os

import pandas as pd

import mitty.benchmarking.dfcols as dfcols

logger = logging.getLogger(__name__)


@click.command()
@click.version_option()
@click.argument('a_fname')
@click.argument('b_fname')
@click.argument('out_prefix')
@click.option('-v', count=True, help='Verbosity level')
def cli(a_fname, b_fname, out_prefix, v):
  """Set operations on two edf or crdf files."""
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  summary = {
    'meta': {
      'A': os.path.basename(a_fname),
      'B': os.path.basename(b_fname)
    }
  }
  for call_type in ['TP', 'FN', 'FP']:
    summary[call_type] = set_operations(
      a_fname, b_fname, call_type=call_type, out_prefix=out_prefix)

  json.dump(summary, open(out_prefix + '-summary.json', 'w'))
  pretty_print_summary(summary, out_prefix + '-summary.txt')


def set_operations(a_fname, b_fname, call_type='FP', out_prefix='setops'):
  """

  :param a_fname:
  :param b_fname:
  :param call_type: {'TP', 'FP', 'FN'}
  :param out_prefix:
  :return:
  """
  sA, f_typeA = get_call_hash_set(a_fname, call_type)
  sB, f_typeB = get_call_hash_set(b_fname, call_type)

  summary = {
    'A': len(sA)
  }

  # A - B
  sA_sB = sA - sB
  summary['A - B'] = len(sA_sB)
  write_out_data(a_fname,
                 dst_fname='{}-a-b-{}.{}.h5'.format(out_prefix, call_type, f_typeA),
                 call_hash_set=sA_sB,
                 f_type=f_typeA)

  # A ^ B
  sAB = sA.intersection(sB)
  summary['A & B'] = len(sAB)
  fname, f_type = (b_fname, f_typeB) if f_typeA == 'edf' and f_typeB == 'crdf' else (a_fname, f_typeA)
  write_out_data(fname,
                 dst_fname='{}-a&b-{}.{}.h5'.format(out_prefix, call_type, f_type),
                 call_hash_set=sAB,
                 f_type=f_type)

  # B - A
  sB_sA = sB - sA
  summary['B - A'] = len(sB_sA)
  write_out_data(b_fname,
                 dst_fname='{}-b-a-{}.{}.h5'.format(out_prefix, call_type, f_typeB),
                 call_hash_set=sB_sA,
                 f_type=f_typeB)

  return summary


def get_call_hash_set(fname, call_type='FP'):
  """

  :param fname:
  :param call_type: {'TP', 'FP', 'FN'}
  :return:
  """
  logger.debug('Loading {} ...'.format(call_type))
  t0 = time.time()

  store = pd.HDFStore(fname, 'r')

  f_type = None
  if 'edf' in store:
    f_type = 'edf'
  elif 'crdf' in store:
    f_type = 'crdf'
  else:
    raise ValueError('Store {} has neither "edf" or "crdf"'.format(fname))

  data = set(
    store.select(
      f_type, columns=['call_hash'],
      where="call_type={}".format(dfcols.call_type[call_type])
    )['call_hash'])

  store.close()

  t1 = time.time()
  logger.debug('... took {:0.5} s'.format(t1 - t0))

  return data, f_type


def write_out_data(src_fname, dst_fname, call_hash_set, f_type):
  logger.debug('Writing out {} ...'.format(dst_fname))
  t0 = time.time()
  pd.read_hdf(src_fname, f_type,
              where=pd.Term('call_hash','=', list(call_hash_set)), mode='r').\
    to_hdf(
    dst_fname,
    f_type, data_columns=dfcols.get_edf_bdf_data_cols().keys(),
    format='t', mode='w', complevel=9, complib='blosc')
  t1 = time.time()
  logger.debug('... took {:0.5} s'.format(t1 - t0))


def pretty_print_summary(summary, fname):
  header = ['TP', 'FN', 'FP']
  with open(fname, 'w') as fp:
    fp.write('A: {}\n'.format(summary['meta']['A']))
    fp.write('B: {}\n'.format(summary['meta']['B']))
    fp.write('------------------------------------\n')
    fp.write('\t'.join([''] + header) + '\n')
    for op in ['A', 'A - B', 'A & B', 'B - A']:
      fp.write('\t'.join([op] + [str(summary[h][op]) for h in header]) + '\n')


if __name__ == '__main__':
  cli()