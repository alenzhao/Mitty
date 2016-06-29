"""Reads in an eval.vcf and turns it into a data frame"""
# http://stackoverflow.com/questions/31090127/pandas-continuously-write-from-function-to-csv
import logging
import gzip
import time
from hashlib import md5


import click
import pandas as pd


logger = logging.getLogger(__name__)


@click.group()
def cli():
  pass


@cli.command('convert')
@click.argument('evalvcf')
@click.argument('outcsv')
@click.option('-v', count=True, help='Verbosity level')
def convert_evcf(evalvcf, outcsv, v):
  """Convert an eval vcf from vcfeval to a dataframe."""
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  read_evcf_into_dataframe(evalvcf).to_csv(outcsv, index=False, compression='gzip' if outcsv.endswith('gz') else None)


@cli.command('compare')
@click.argument('file_a')
@click.argument('file_b')
@click.argument('file_a_and_b')
@click.argument('file_a_only')
@click.argument('file_b_only')
def set_operations_evcf(file_a, file_b, file_a_and_b, file_a_only, file_b_only):
  """Perform intersection and difference operations on the variant calls."""
  dfA = pd.read_csv(file_a, compression='gzip' if file_a.endswith('gz') else None).set_index(['chrom', 'pos', 'ref', 'alt', 'truth', 'query'])
  dfB = pd.read_csv(file_a, compression='gzip' if file_a.endswith('gz') else None).set_index(['chrom', 'pos', 'ref', 'alt', 'truth', 'query'])

  # Intersect
  dfA[dfA.index.isin(dfB.index)].to_csv(file_a_and_b, compression='gzip' if file_a_and_b.endswith('gz') else None)

  # A - B
  dfA[~dfA.index.isin(dfB.index)].to_csv(file_a_only, compression='gzip' if file_a_only.endswith('gz') else None)

  # B - A
  dfB[~dfB.index.isin(dfA.index)].to_csv(file_b_only, compression='gzip' if file_b_only.endswith('gz') else None)


def get_header_row_count_for_vcf(fname):
  with gzip.open(fname) as fp:
    for cnt, line in enumerate(fp):
      if not line.startswith('##'):
        return cnt
  return None


def get_contig_dict(fname):
  seq_dict = {}
  seq_cntr = 1
  with gzip.open(fname) as fp:
    for line in fp:
      if line.startswith('##contig=<ID'):
        x = line.split(',', 1)
        seq_dict[x[0][13:]] = seq_cntr
        seq_cntr += 1
      if not line.startswith('##'):
        break
  return seq_dict


# This should be in a central library
def call_hash(row):
  return int(md5(','.join([row['#CHROM'], str(row['POS']), row['REF'], row['ALT']])).hexdigest()[:8], 16)


def call_type_encoding(row):
  if row['truth'] == 'TP':
    return 0
  elif row['truth'] == 'FN':
    return 1
  else:
    return 2  # FP


def read_evcf_into_dataframe(fname):

  t0 = time.time()
  hdr_cnt = get_header_row_count_for_vcf(fname)
  seq_dict = get_contig_dict(fname)
  evcf = pd.read_csv(fname, dtype={'#CHROM': 'S10'}, skiprows=hdr_cnt, sep='\t', compression='gzip', engine='c')
  t1 = time.time()
  logger.debug('Loaded eval csv into basic data frame ({:0.3} s) ({} rows)'.format(t1 - t0, len(evcf)))

  df = pd.DataFrame()
  df['call_chrom'] = evcf['#CHROM']
  df['call_pos'] = evcf['POS']
  df['ref'] = evcf['REF']
  df['alt'] = evcf['ALT']
  t1_5 = time.time()
  logger.debug('Copied selected columns ({:0.3} s)'.format((t1_5 - t1)))

  df['call_hash'] = evcf.apply(call_hash, axis=1)
  t2 = time.time()
  logger.debug('Computed hash for call ({:0.3} s)'.format((t2 - t1_5)))

  df['ROC_thresh'], df['truth'], df['truthGT'], df['query'], df['queryGT'] = parse_call_column(evcf)
  t3 = time.time()
  logger.debug('Parsed call columns using format string ({:0.3} s)'.format((t3 - t2)))

  df['call_type'] = df.apply(call_type_encoding, axis=1)
  t3_5 = time.time()
  logger.debug('Computed call type ({:0.3} s)'.format((t3_5 - t3)))

  df['variant_size'] = evcf.apply(variant_size, axis=1)
  t4 = time.time()
  logger.debug('Computed variant size ({:0.3} s)'.format((t4 - t3_5)))

  # df['pos_stop'] = df['pos'] - df['variant_size'].apply(lambda x: min(x, 0))
  evcf['pos_stop'] = df['call_pos'] - df['variant_size'].apply(lambda x: min(x, 0))
  evcf['int_chrom'] = evcf.apply(lambda row: seq_dict[row['#CHROM']], axis=1)
  t5 = time.time()
  logger.debug('Converted alphanumeric chrom code to integer seq_id ({:0.3} s)'.format((t5 - t4)))

  df['call_p1'] = evcf.apply(lambda row: (row['int_chrom'] << 29) | row['POS'], axis=1)
  df['call_p2'] = evcf.apply(lambda row: (row['int_chrom'] << 29) | row['pos_stop'], axis=1)
  t6 = time.time()
  logger.debug('Computed compressed coordinates ({:0.3} s)'.format((t6 - t5)))

  logger.debug('Took {:0.3} s to process eval.vcf'.format(t6 - t0))

  return df.sort_values('call_p1')


# TP:
#   #CHROM                              10
#   POS                            6232682
#   ID                                   .
#   REF                                  C
#   ALT                                  A
#   QUAL                                 .
#   FILTER                               .
#   INFO                        BS=6232682
#   FORMAT          GT:BD:BK:QQ:BI:BVT:BLT
#   TRUTH           0|1:TP:gm:.:tv:SNP:het
#   QUERY     0/1:TP:gm:165.129:tv:SNP:het
#
#
# FP:
#   #CHROM                                 1
#   POS                              1012309
#   ID                                     .
#   REF                                    T
#   ALT                                    C
#   QUAL                                   .
#   FILTER                                 .
#   INFO                                BS=1
#   FORMAT            GT:QQ:BD:BK:BI:BVT:BLT
#   TRUTH            .:.:.:.:.:NOCALL:nocall
#   QUERY     1/1:28.4947:FP:.:ti:SNP:homalt
#
#
# FN:
#   #CHROM                                                   10
#   POS                                                 6228641
#   ID                                                        .
#   REF       AGGAAAACTTACTCAGAGCCTAATCTGTAAACCAAACCTATGGGGA...
#   ALT                                                       A
#   QUAL                                                      .
#   FILTER                                                    .
#   INFO                                             BS=6228642
#   FORMAT                                  GT:BD:BK:BI:BVT:BLT
#   TRUTH                           0|1:FN:.:d16_plus:INDEL:het
#   QUERY                                 .:.:.:.:NOCALL:nocall


# column splitting from http://www.swegler.com/becky/blog/2014/08/06/useful-pandas-snippets/
def parse_call_column(evcf):
  """Returns multiple columns:
  ROC_thresh, truth, truthGT, query, queryGT

  """
  def parse(row):
    keys = row['FORMAT'].split(':')
    tv = {k: v for k, v in zip(keys, row['TRUTH'].split(':'))}
    qv = {k: v for k, v in zip(keys, row['QUERY'].split(':'))}
    return qv.get('QQ', None), tv['BD'], tv['GT'], qv['BD'], qv['GT']

  return zip(*evcf.apply(parse, axis=1))


def variant_size(val):
  return len(val['ALT']) - len(val['REF'])


if __name__ == '__main__':
  cli()