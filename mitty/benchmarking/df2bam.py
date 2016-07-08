import click
import logging
import time
from hashlib import md5
import os
from copy import deepcopy

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

  df = load_df(df_fname)
  qnames = set(df['qname_hash'])
  writeout_bam(qnames, df, origbam, outbam)
  # verify_bam(df_fname, outbam)


def load_df(df_fname):
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

  df = store.select(f_type)
  # from IPython import embed; embed()
  # exit()
  # data = set(store.select_column(f_type, 'qname_hash'))

  store.close()

  t1 = time.time()
  logger.debug('... took {:0.5} s'.format(t1 - t0))

  return df


def writeout_bam(qnames, df, origbam, outbam):
  logger.debug('Writing out BAM ... ')
  t0 = time.time()
  fp_in = pysam.AlignmentFile(origbam, 'rb')

  outbam_per = os.path.splitext(outbam)[0] + '.per.bam'
  fp_out = pysam.AlignmentFile(outbam, 'wb', header=fp_in.header)
  fp_out_per = pysam.AlignmentFile(outbam_per, 'wb', header=fp_in.header)

  notification_interval = int(0.1 * len(qnames))
  ctr = 0
  for r in fp_in:
    q_h = int(md5(r.qname).hexdigest()[:8], 16)
    if q_h in qnames:
      row = df[df['qname_hash'] == q_h].iloc[0]
      mate = 'm1' if r.is_read1 else 'm2'  # Non paired end?
      Xd = row[mate + '_d_error']
      r.set_tag('Xd', Xd)
      fp_out.write(r)  # This writes out the original read with only the Xd tag added

      p1 = row[mate + '_correct_p1']
      c_chrom = p1 >> 29
      c_pos = p1 & 0x1fffffff
      r.reference_id = c_chrom - 1
      r.pos = c_pos
      r.is_unmapped = not (row[mate + '_mapped'] & 0b10) >> 1
      fp_out_per.write(r)  # This writes out the read placed in the correct location

      ctr += 1
      if ctr % notification_interval == 0:
        t1 = time.time()
        logger.debug('Wrote out {} reads {:0.5} sec'.format(ctr, t1 - t0))

  fp_in.close()
  fp_out.close()
  fp_out_per.close()

  t1 = time.time()
  logger.debug('... took {:0.5} s'.format(t1 - t0))

  t0 = time.time()
  mio.sort_and_index_bam(outbam)
  mio.sort_and_index_bam(outbam_per)
  t1 = time.time()
  logger.debug('Sort and indexed BAM fragments in {:2.2f}s'.format(t1 - t0))


# def verify_bam(df_fname, outbam):
#   print('Verifying BAM')
#   store = pd.HDFStore(df_fname, 'r')
#   my_qnames = set(store.select('crdf', where='call_hash=2730745348')['qname_hash'])
#   df = store.select('crdf', where='call_hash=2730745348')
#
#   print('{} qnames'.format(len(my_qnames)))
#
#   fp_in = pysam.AlignmentFile(outbam, 'rb')
#   for r in fp_in:
#     q_h = int(md5(r.qname).hexdigest()[:8], 16)
#     if q_h in my_qnames:
#       print(r.qname, r.pos, bin(df[df['qname_hash'] == q_h]['align_flag'].values[0]), df[df['qname_hash'] == q_h][['m1_d_error', 'm2_d_error']].values[0])
#       # my_qnames.remove(q_h)
#       # print(len(my_qnames))
#       # if len(my_qnames) == 0: break


if __name__ == '__main__':
  cli()