import click
import logging
import time
import gzip

import pysam  # For tabix/bgzip compression
import pandas as pd

import mitty.lib.mio as mio  # For the bam sort and index function
from mitty.benchmarking.evalvcf2df import get_header_row_count_for_vcf, get_contig_dict
import mitty.benchmarking.dfcols as dfcols


logger = logging.getLogger(__name__)


@click.command()
@click.version_option()
@click.argument('df_fname')
@click.argument('origvcf')
@click.argument('outvcf')
@click.option('-v', count=True, help='Verbosity level')
def cli(df_fname, origvcf, outvcf, v):
  """Construct a BAM fragment based on a bf or crdf file. Be careful of memory
  """
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  call_hashes = get_call_hash_set(df_fname)
  writeout_vcf(call_hashes, origvcf, outvcf)


def get_call_hash_set(fname):
  """

  :param fname:
  :return:
  """
  logger.debug('Loading df ...')
  t0 = time.time()

  store = pd.HDFStore(fname, 'r')

  f_type = None
  if 'edf' in store:
    f_type = 'edf'
  elif 'crdf' in store:
    f_type = 'crdf'
  else:
    raise ValueError('Store {} has neither "edf" or "crdf"'.format(fname))

  data = set(store.select_column(f_type, 'call_hash'))

  store.close()

  t1 = time.time()
  logger.debug('... took {:0.5} s'.format(t1 - t0))

  return data


# For somereason, gzipped
def writeout_vcf(call_hashes, origvcf, outvcf):
  t0 = time.time()

  evcf = load_evcf(origvcf)
  t1 = time.time()
  logger.debug('Loaded eval csv into basic data frame ({:0.3} s) ({} rows)'.format(t1 - t0, len(evcf)))

  evcf_call_hashes = evcf.apply(dfcols.call_hash, axis=1)
  t2 = time.time()
  logger.debug('Computed hash for call ({:0.3} s)'.format((t2 - t1)))

  idx = evcf_call_hashes.isin(call_hashes)

  logger.debug('Writing out VCF ...')

  bgzip_me, unzipped_file = (True, outvcf[:-3]) if outvcf.endswith('.gz') else (False, outvcf)

  # Write out header
  with gzip.open(origvcf, 'r') as fp_in, open(unzipped_file, 'w') as fp_out:
    for cnt, line in enumerate(fp_in):
      if line.startswith('##'):
        fp_out.write(line)
      else:
        break

  evcf[idx].to_csv(unzipped_file, mode='a', index=False, sep='\t')

  if bgzip_me:
    pysam.tabix_compress(unzipped_file, outvcf, force=True)

  t3 = time.time()
  logger.debug('... took {:0.5} s'.format(t3 - t2))


def load_evcf(origvcf):
  hdr_cnt = get_header_row_count_for_vcf(origvcf)
  seq_dict = get_contig_dict(origvcf)
  chrom_id_max_len = max([len(v) for v in seq_dict.keys()])
  chrom_id_fmt = 'a{}'.format(chrom_id_max_len)
  return pd.read_csv(origvcf, dtype={'#CHROM': chrom_id_fmt}, skiprows=hdr_cnt, sep='\t', compression='gzip', engine='c')


if __name__ == '__main__':
  cli()