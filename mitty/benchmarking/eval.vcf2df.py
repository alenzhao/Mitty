"""Reads in an eval.vcf and turns it into a data frame"""
# http://stackoverflow.com/questions/31090127/pandas-continuously-write-from-function-to-csv

import gzip

import click
import pandas as pd


@click.command()
@click.argument('evalvcf')
@click.argument('outcsv')
# @click.option('--block-size', default=5000)
def cli(evalvcf, outcsv):
  parse_evcf(read_evcf_into_dataframe(evalvcf)).to_csv(outcsv, index=False, compression='gzip' if outcsv.endswith('gz') else None)



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


def parse_call_column(val):
  """TP, FP, FN"""
  if val['TRUTH'][0] == '.':
    return 'FP'
  elif val['QUERY'][0] == '.':
    return 'FN'
  else:
    return 'TP'


def variant_size(val):
  return len(val['ALT']) - len(val['REF'])


def parse_evcf(evcf):
  df = pd.DataFrame()
  df['chrom'] = evcf['#CHROM']
  df['pos'] = evcf['POS']
  df['ref'] = evcf['REF']
  df['alt'] = evcf['ALT']
  df['call'] = evcf.apply(parse_call_column, axis=1)
  df['variant_category'] = evcf.apply(variant_size, axis=1)
  return df

if __name__ == '__main__':
  cli()