"""Reads in an eval.vcf and turns it into a data frame"""
# http://stackoverflow.com/questions/31090127/pandas-continuously-write-from-function-to-csv

import gzip

import click
import pandas as pd


@click.group()
def cli():
  pass
#  parse_evcf(read_evcf_into_dataframe(evalvcf)).to_csv(outcsv, index=False, compression='gzip' if outcsv.endswith('gz') else None)


@cli.command('convert')
@click.argument('evalvcf')
@click.argument('outcsv')
def convert_evcf(evalvcf, outcsv):
  """Convert an eval vcf from vcfeval to a dataframe."""
  parse_evcf(read_evcf_into_dataframe(evalvcf)).to_csv(outcsv, index=False, compression='gzip' if outcsv.endswith('gz') else None)


def set_operations_evcf(csvA, csvB):
  pass





  # from IPython import embed; embed()
  # read_evcf_into_dataframe(evalvcf).to_csv(out_csv, index=False)

  # for cols in fastq_reader_ip(fname, block_size):
  #   df = pd.DataFrame(cols,
  #                     columns=['serial', 'qname', 'chrom', 'phase', 'strand', 'pos', 'rlen', 'cigar'])
  #   df.to_csv(out_csv, mode='a', index=False)


def get_header_row_count_for_vcf(fname):
  with gzip.open(fname) as fp:
    for cnt, line in enumerate(fp):
      if not line.startswith('##'):
        return cnt
  return None


def read_evcf_into_dataframe(fname):
  hdr_cnt = get_header_row_count_for_vcf(fname)
  evcf = pd.read_csv(fname, skiprows=hdr_cnt, sep='\t', compression='gzip', engine='c')
  return evcf


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
  #_fmt, _truth, _query = evcf['FORMAT'], evcf['TRUTH'], evcf['QUERY']

  def parse(row):
    keys = row['FORMAT'].split(':')
    tv = {k: v for k, v in zip(keys, row['TRUTH'].split(':'))}
    qv = {k: v for k, v in zip(keys, row['QUERY'].split(':'))}
    return qv.get('QQ', None), tv['BD'], tv['GT'], qv['BD'], tv['GT']

  return zip(*evcf.apply(parse, axis=1))


def variant_size(val):
  return len(val['ALT']) - len(val['REF'])


def parse_evcf(evcf):
  # from IPython import embed; embed()
  # exit()

  df = pd.DataFrame()
  df['chrom'] = evcf['#CHROM']
  df['pos'] = evcf['POS']
  df['ref'] = evcf['REF']
  df['alt'] = evcf['ALT']
  df['ROC_thresh'], df['truth'], df['truthGT'], df['query'], df['queryGT'] = parse_call_column(evcf)
  df['variant_category'] = evcf.apply(variant_size, axis=1)
  return df.set_index(['chrom', 'pos', 'ref', 'alt'])

if __name__ == '__main__':
  cli()