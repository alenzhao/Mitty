"""Read in an aligner (qname sorted) BAM and write out a CSV with pairs of reads with the aligned position,
correct position and tags.

Continuous writing from http://stackoverflow.com/questions/31090127/pandas-continuously-write-from-function-to-csv
"""
import logging
import time
from multiprocessing import Pool
import itertools

import click
import pysam
import pandas as pd

logger = logging.getLogger(__name__)

read_info = [
  'd_error',   # Alignment error metric
  'a_mapped',  # Aligned mapped status
  'a_pos',     # Aligned pos = chrom << 29 | pos
  'c_mapped',  # Correct mapped status
  'c_pos',     # Correct pos
  'MQ'
]

gral_tags = [
  'XG',
  'YG',
  'Yg',
  'XI',
  'XB',
  'XE',
  'YX',
  'Yx',
  'YS',
  'YA',
  'YQ',
  'Ym',
  'UQ'
]


@click.command()
@click.argument('qnamesortedbam')
@click.argument('outcsv')
@click.option('--paired-reads/--single-end-reads', default=True, help='Are these paired end or single end reads')
@click.option('--simulated-reads/--real-reads', default=True, help='Are these simulated from Mitty? If so, parse qnames for truth positions')
@click.option('-v', count=True, help='Verbosity level')
def cli(qnamesortedbam, outcsv, paired_reads, simulated_reads, v):
  """Associate reads in a BAM with variant calls."""
  # print(alignerbam)
  # print(evalvcfdf)
  # print(outcsv)
  # print(paired_reads)
  # print(simulated_reads)

  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  logger.warning('This code only works for paired end simulated files ...')
  logger.warning('It is trvial to convert it to work for data with no qname information and SE')

  process_bam(qnamesortedbam, outcsv)



# for cnt, r1 in enumerate(fp):
#   r2 = next(fp)
#
#    assert r1.qname == r2.qname
#    r_d = [parse_qname(r1), parse_qname(r2)]
#    df_rows.append([r1.qname] + [str(r_dd.get(k, '')) for r_dd in r_d for k in read_info + gral_tags])
##   out_fp.write(','.join([r1.qname] + [str(r_dd.get(k, '')) for r_dd in r_d for k in read_info + gral_tags]) + '\n')
  # for r1r2 in get_read_pairs(fp):
  #   df_rows.append(parse_pair(r1r2))


def process_bam(qnamesortedbam, out_csv, update_every_n_templates=10000):
  out_fp = open(out_csv, 'w')

  columns = ['qname'] + [m + t for m in ['m1_', 'm2_'] for t in read_info + gral_tags]
  out_fp.write(','.join(columns) + '\n')

  fp = pysam.AlignmentFile(qnamesortedbam)

  pool = Pool(4)
  block_size = 1000
  pp_fp = get_read_pairs(fp)
  while 1:
    df_rows = pool.map(parse_pair, itertools.islice(pp_fp, block_size))
    pd.DataFrame(df_rows, columns=columns).to_csv(out_csv, index=False, mode='a', compression='gzip' if out_csv.endswith('gz') else None)


def get_read_pairs(fp, update_every_n_templates=10000):
  t0 = time.time()
  for cnt, r1 in enumerate(fp):
    r2 = next(fp)

    if cnt % update_every_n_templates == 0:
      t1 = time.time()
      logger.debug('{} templates processed ({} templates/sec)'.format(cnt, cnt/(t1 - t0)))
      if cnt > 30000:
        break

    yield (r1, r2)


def parse_pair(r1r2):
  return [r1r2[0].qname] + [str(r_dd.get(k, '')) for r_dd in [parse_qname(r1r2[0]), parse_qname(r1r2[1])] for k in read_info + gral_tags]


def parse_qname(read):

  # Took out unpaired read handling for speed

  # if read.is_paired:
  #   if read.is_read1:
  #     rs, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m = read.qname.split('|')
  #   else:
  #     rs, chrom, cpy, ro_m, pos_m, rl_m, cigar_m, ro, pos, rl, cigar = read.qname.split('|')
  # else:
  #   rs, chrom, cpy, ro, pos, rl, cigar = read.qname.split('|')[:7]

  # Parse qname
  if read.is_read1:
    rs, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m = read.qname.split('|')
  else:
    rs, chrom, cpy, ro_m, pos_m, rl_m, cigar_m, ro, pos, rl, cigar = read.qname.split('|')

  chrom, pos = int(chrom), int(pos)

  a_cigar = read.cigarstring

  # WARNING: this destructively changes the orignal cigar string of the read, so we preseve it above
  # We do this to reuse the AlignedSegment object, which is expensive to create
  read.cigarstring = cigar
  cigar_ops = read.cigartuples

  c_mapped = 1
  # Parse alignment
  long_insert = 0
  if read.is_unmapped:
    a_mapped = 0
    d = 1000000000
  else:
    a_mapped = 1
    if read.reference_id != chrom - 1:
      d = 1000000000
    else:  # Analyze the correctness by checking each breakpoint
      # Corner case, our special cigar for indicating reads inside an insertion
      # We use S or I for this
      if cigar_ops[0][0] in [1, 4] and len(cigar_ops) == 1:  # S,I
        d = read.pos - pos
        c_mapped = 0  # In a long insert
      else:  # Go through breakpoints
        correct_pos = pos
        d = read.pos - correct_pos
        for op, cnt in cigar_ops:
          if op in [0, 7, 8]:  # M, =, X
            this_d = read.pos - correct_pos
            if abs(this_d) < abs(d):
              d = this_d
            correct_pos += cnt
          elif op == 2:  # D
            correct_pos += cnt

  a_pos = ((read.reference_id + 1) << 29) | read.pos
  c_pos = (chrom << 29) | pos

  r_dict = {
    'd_error': d,  # Fill this out with proper metric
    'qname': read.qname,
    'mate': 0 if read.is_read1 else 1,
    'a_mapped': a_mapped,
    'a_pos': a_pos,
    'a_cigar': a_cigar,
    'c_mapped': c_mapped,
    'c_pos': c_pos,
    'c_cigar': cigar,
    'MQ': read.mapping_quality,
  }
  for k, v in read.get_tags():
    r_dict[k] = v
  return r_dict


if __name__ == '__main__':
  cli()