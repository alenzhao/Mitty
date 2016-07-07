import click
import logging
import time
from hashlib import md5

import pandas as pd
import pysam

import mitty.lib.mio as mio  # For the bam sort and index function

logger = logging.getLogger(__name__)


@click.command()
@click.version_option()
@click.argument('df_fname')
@click.argument('origbam')
@click.argument('outbam')
@click.option('-v', count=True, help='Verbosity level')
def cli(df_fname, origbam, outbam, v):
  """Construct a BAM fragment based on a bf or crdf file. Be careful of memory
  """
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  qnames = get_qnames(df_fname)
  writeout_bam(qnames, origbam, outbam)


def get_qnames(df_fname):
  logger.debug('Loading df ...')
  t0 = time.time()

  store = pd.HDFStore(df_fname, 'r')

  f_type = None
  if 'bdf' in store:
    f_type = 'bdf'
  elif 'crdf' in store:
    f_type = 'crdf'
  else:
    raise ValueError('Store {} has neither "bdf" or "crdf"'.format(df_fname))

  data = set(store.select_column(f_type, 'qname_hash'))

  store.close()

  t1 = time.time()
  logger.debug('... took {:0.5} s'.format(t1 - t0))

  return data


def writeout_bam(qnames, origbam, outbam):
  logger.debug('Writing out BAM ... ')
  t0 = time.time()
  fp_in = pysam.AlignmentFile(origbam, 'rb')
  fp_out = pysam.AlignmentFile(outbam, 'wb', header=fp_in.header)
  notification_interval = int(0.1 * len(qnames))
  ctr = 0
  for r in fp_in:
    if int(md5(r.qname).hexdigest()[:8], 16) in qnames:
      fp_out.write(r)
      ctr += 1
      if ctr % notification_interval == 0:
        t1 = time.time()
        logger.debug('Wrote out {} reads {:0.5} sec'.format(ctr, t1 - t0))

  fp_in.close()
  fp_out.close()

  t1 = time.time()
  logger.debug('... took {:0.5} s'.format(t1 - t0))

  t0 = time.time()
  mio.sort_and_index_bam(outbam)
  t1 = time.time()
  logger.debug('Sort and indexed bad BAM in {:2.2f}s'.format(t1 - t0))


if __name__ == '__main__':
  cli()