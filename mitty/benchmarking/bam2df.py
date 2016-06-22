"""Read in an aligner (qname sorted) BAM and write out a CSV with pairs of reads with the aligned position,
correct position and tags.

Continuous writing from http://stackoverflow.com/questions/31090127/pandas-continuously-write-from-function-to-csv


python ~/Code/Mitty/mitty/benchmarking/bam2df.py qname-run00001.S.pe.100x500.perfect.g-G1.t-0.8.17.bam bam.csv -v


"""
import logging
import time
from multiprocessing import Pool
import subprocess

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
@click.option('--block-size', default=100000, help='Block size')
@click.option('-t', default=4, help='Threads')
@click.option('-v', count=True, help='Verbosity level')
def cli(qnamesortedbam, outcsv, paired_reads, simulated_reads, block_size, t, v):
  """Associate reads in a BAM with variant calls."""
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  logger.warning('This code only works for paired end simulated files ...')
  logger.warning('It is trvial to convert it to work for data with no qname information and SE')

  process_bam_parallel(qnamesortedbam, outcsv, block_size, threads=t)


def process_bam_parallel(qnamesortedbam, out_csv, block_size=10000, threads=4):
  """Header is in 'tmp_hdr.csv' partials are in tmp_1.csv, tmp_2.csv ...."""
  prefix = 'temp_bam2df'

  columns = ['qname'] + [m + t for m in ['m1_', 'm2_'] for t in read_info + gral_tags]

  # Write header
  hdr_fname = prefix + '_hdr.csv'
  with open(hdr_fname, 'w') as out_fp:
    out_fp.write(','.join(columns) + '\n')

  offsets = break_bam(qnamesortedbam, block_size)
  fnames = [prefix + '_{}.csv'.format(n) for n in range(len(offsets))]

  p = Pool(threads)
  t_total = 0
  t0 = time.time()
  for fn in p.imap_unordered(process_bam_section, [(qnamesortedbam, off, fn, block_size) for off, fn in zip(offsets, fnames)]):
    print fn
    t1 = time.time()
    t_total += block_size
    logger.debug('{} templates processed ({} templates/sec)'.format(t_total, t_total/(t1 - t0)))

  # This will only work on unix like systems
  with open(out_csv, 'w') as fp:
    subprocess.call(['cat', hdr_fname] + fnames, stdout=fp)

  subprocess.call('rm {}*.csv'.format(prefix), shell=True)


def break_bam(qnamesortedbam, block_size):
  """Go through the BAM and mark-off the offsets,

  :param qnamesortedbam:
  :param block_size:
  :return:
  """
  t_total = 0
  t0 = time.time()
  offsets = [off for off in get_bam_sections(qnamesortedbam, block_size)]
  t1 = time.time()
  logger.debug('Pass 1: Find points to split BAM {} sec to run through BAM'.format(t1 - t0))

  return offsets


def get_bam_sections(qnamesortedbam, block_size):
  fp = pysam.AlignmentFile(qnamesortedbam)
  for cnt, _ in enumerate(fp):
    next(fp)
    if cnt % block_size == 0:
      yield fp.tell()


# args = qnamesortedbam, template_start, template_end
def process_bam_section(args):
  qnamesortedbam, offset, fname, block_size = args
  fp = pysam.AlignmentFile(qnamesortedbam)
  fp.seek(offset)

  out_fp = open(fname, 'w')
  for row in [parse_pair(r1r2) for r1r2 in get_read_pairs(fp, block_size)]:
    out_fp.write(','.join(row) + '\n')

  return fname

def get_read_pairs(fp, block_size):
  for cnt, r1 in enumerate(fp):
    if cnt == block_size:
      break
    r2 = next(fp)
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




def process_bam(qnamesortedbam, out_csv, block_size=10000):
  out_fp = open(out_csv, 'w')

  columns = ['qname'] + [m + t for m in ['m1_', 'm2_'] for t in read_info + gral_tags]
  out_fp.write(','.join(columns) + '\n')

  t_total = 0
  t0 = time.time()
  for offset in get_bam_sections(qnamesortedbam, block_size):
    df_rows = process_bam_section((qnamesortedbam, offset, block_size))
    pd.DataFrame(df_rows, columns=columns).to_csv(out_csv, index=False, mode='a', compression='gzip' if out_csv.endswith('gz') else None)
    t1 = time.time()
    t_total += block_size
    logger.debug('{} templates processed ({} templates/sec)'.format(t_total, t_total/(t1 - t0)))


# https://groups.google.com/forum/#!topic/pysam-user-group/bBdqn7DkVtE
# Using the multiprocessing module should work, but to be on the safe
# side, open a separate file for read access in each worker
# process. That is, instead of passing an instance of Samfile to a
# subprocess, pass the filename instead and create a new samfile =
# Samfile( ) in the function that is run in parallel.
#
# There are definitely side-effects when multi-threading and
# possibly ones when multi-processing. Generally it is best to avoid
# accessing a bam-file through the same instance of a Samfile object
# from multiple threads/processes.

def process_bam_parallel2(qnamesortedbam, out_csv, block_size=10000):
  out_fp = open(out_csv, 'w')

  columns = ['qname'] + [m + t for m in ['m1_', 'm2_'] for t in read_info + gral_tags]
  out_fp.write(','.join(columns) + '\n')

  t_total = 0
  p = Pool(8)
  t0 = time.time()
  for df_rows in p.imap(
    process_bam_section,
    ((qnamesortedbam, offset, block_size) for offset in get_bam_sections(qnamesortedbam, block_size))):
    pd.DataFrame(df_rows, columns=columns).to_csv(out_csv, index=False, mode='a', compression='gzip' if out_csv.endswith('gz') else None)
    t1 = time.time()
    t_total += block_size
    logger.debug('{} templates processed ({} templates/sec)'.format(t_total, t_total/(t1 - t0)))


def process_bam_parallel3(qnamesortedbam, out_csv, block_size=100000):
  # out_fp = open(out_csv, 'w')
  #
  columns = ['qname'] + [m + t for m in ['m1_', 'm2_'] for t in read_info + gral_tags]
  # out_fp.write(','.join(columns) + '\n')

  t_total = 0
  # p = Pool(2)
  t0 = time.time()
  offsets = [off for off in get_bam_sections(qnamesortedbam, block_size)]
  t1 = time.time()
  logger.debug('Pass 1: Find points to split BAM {} sec to run through BAM'.format(t1 - t0))


  # for df_rows in p.imap(
  #   process_bam_section,
  #   ((qnamesortedbam, offset, block_size) for offset in get_bam_sections(qnamesortedbam, block_size))):
  #   pd.DataFrame(df_rows, columns=columns).to_csv(out_csv, index=False, mode='a', compression='gzip' if out_csv.endswith('gz') else None)
  #   t1 = time.time()
  #   t_total += block_size
  #   logger.debug('{} templates processed ({} templates/sec)'.format(t_total, t_total/(t1 - t0)))

  # for df_rows in itertools.imap(
  #   process_bam_section,
  #   ((qnamesortedbam, offset, block_size) for offset in get_bam_sections(qnamesortedbam, block_size))):
  #   pd.DataFrame(df_rows, columns=columns).to_csv(out_csv, index=False, mode='a', compression='gzip' if out_csv.endswith('gz') else None)
  #   t1 = time.time()
  #   t_total += block_size
  #   logger.debug('{} templates processed ({} templates/sec)'.format(t_total, t_total/(t1 - t0)))


def process_bam_parallel4(qnamesortedbam, out_csv, block_size=1000):
  out_fp = open(out_csv, 'w')

  columns = ['qname'] + [m + t for m in ['m1_', 'm2_'] for t in read_info + gral_tags]
  out_fp.write(','.join(columns) + '\n')

  t_total = 0
  p = Pool()

  fp = pysam.AlignmentFile(qnamesortedbam)

  t0 = time.time()
  for rp_l in get_read_pairs_in_blocks(fp, block_size):
    df_rows = p.map(parse_pair, rp_l, chunksize=block_size/4)
    # df_rows = map(parse_pair, rp_l)
    pd.DataFrame(df_rows, columns=columns).to_csv(out_csv, index=False, mode='a', compression='gzip' if out_csv.endswith('gz') else None)
    t1 = time.time()
    t_total += block_size
    logger.debug('{} templates processed ({} templates/sec)'.format(t_total, t_total/(t1 - t0)))

  # for df_rows in p.imap_unordered(process_pairs, (rpp for rpp in get_read_pairs_in_blocks(fp, block_size))):
  #   pd.DataFrame(df_rows, columns=columns).to_csv(out_csv, index=False, mode='a', compression='gzip' if out_csv.endswith('gz') else None)
  #   t1 = time.time()
  #   t_total += block_size
  #   logger.debug('{} templates processed ({} templates/sec)'.format(t_total, t_total/(t1 - t0)))



def get_read_pairs_in_blocks(fp, block_size):
  rpp = []
  for cnt, r1 in enumerate(fp):
    r2 = next(fp)
    rpp += [(r1, r2)]
    if (cnt + 1) % block_size == 0:
      yield rpp
      rpp = []

  yield rpp


def get_all_read_pairs(fp):
  for r1 in fp:
    r2 = next(fp)
    yield tuple([r1, r2])


def process_pairs(r1r2_g):
  return [parse_pair(r1r2) for r1r2 in r1r2_g]


if __name__ == '__main__':
  cli()