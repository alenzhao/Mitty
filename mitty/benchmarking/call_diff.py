import click
import logging
import time

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
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  for call_type in ['TP', 'FN', 'FP']:
    set_operations(a_fname, b_fname, call_type=call_type, out_prefix=out_prefix)


def set_operations(a_fname, b_fname, call_type='FP', out_prefix='setops'):
  """

  :param a_fname:
  :param b_fname:
  :param call_type: {'TP', 'FP', 'FN'}
  :param out_prefix:
  :return:
  """
  sA, f_typeA = get_qname_set(a_fname, call_type)
  sB, f_typeB = get_qname_set(b_fname, call_type)

  # A - B
  write_out_data(a_fname,
                 dst_fname='{}-a-b-{}.{}.h5'.format(out_prefix, call_type, f_typeA),
                 call_hash_set=sA - sB,
                 f_type=f_typeA)

  # A ^ B
  fname, f_type = (b_fname, f_typeB) if f_typeA == 'edf' and f_typeB == 'crdf' else (a_fname, f_typeA)
  # if f_typeA == 'edf' and f_typeB == 'crdf':
  #   fname = b_fname
  #   f_type = f_typeB
  # else:
  #   fname = a_fname
  #   f_type = f_typeA
  write_out_data(fname,
                 dst_fname='{}-a&b-{}.{}.h5'.format(out_prefix, call_type, f_type),
                 call_hash_set=sA.intersection(sB),
                 f_type=f_type)

  # B - A
  write_out_data(b_fname,
                 dst_fname='{}-b-a-{}.{}.h5'.format(out_prefix, call_type, f_typeB),
                 call_hash_set=sB - sA,
                 f_type=f_typeB)


def get_qname_set(fname, call_type='FP'):
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
    )['call_hash'].unique())

  store.close()

  t1 = time.time()
  logger.debug('... took {:10.2} s'.format(t1 - t0))

  return data, f_type


def write_out_data(src_fname, dst_fname, call_hash_set, f_type):

  logger.debug('Writing out {} ...'.format(dst_fname))
  t0 = time.time()
  pd.read_hdf(src_fname, f_type,
              where=pd.Term('call_hash','=', list(call_hash_set)), mode='r').\
    to_hdf(dst_fname, f_type, format='t', mode='w', complevel=9, complib='blosc')
  t1 = time.time()
  logger.debug('... took {:10.2} s'.format(t1 - t0))


if __name__ == '__main__':
  cli()