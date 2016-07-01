"""Read in an aligner (qname sorted) BAM and write out a CSV with pairs of reads with the aligned position,
correct position and tags.


python ~/Code/Mitty/mitty/benchmarking/bam2df.py qname-run00001.S.pe.100x500.perfect.g-G1.t-0.8.17.bam read-call.csv -v -t 2

"""
import logging
import time
from multiprocessing import Pool
import os
import subprocess
from hashlib import md5

import numpy as np
import pandas as pd
import click
import pysam

import mitty.benchmarking.dfcols as dfcols

logger = logging.getLogger(__name__)


@click.command()
@click.argument('qnamesortedbam')
@click.argument('outh5')
@click.option('--paired-reads/--single-end-reads', default=True, help='Are these paired end or single end reads')
@click.option('--simulated-reads/--real-reads', default=True, help='Are these simulated from Mitty? If so, parse qnames for truth positions')
@click.option('--block-size', default=100000, help='Block size')
@click.option('--max-templates', type=int, help='Max templates to process (for debugging)')
@click.option('-t', default=4, help='Threads')
@click.option('-v', count=True, help='Verbosity level')
def cli(qnamesortedbam, outh5, paired_reads, simulated_reads, block_size, max_templates, t, v):
  """Associate reads in a BAM with variant calls."""
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  logger.warning('This code only works for paired end simulated files ...')
  logger.warning('It is trvial to convert it to work for data with no qname information and SE')

  if os.path.exists(outh5):
    logger.warning('Removing existing file {}'.format(outh5))
    os.remove(outh5)

  if max_templates is not None:
    logger.debug('Will process only {:d} templates'.format(max_templates))

  process_bam_parallel(qnamesortedbam, outh5, block_size, threads=t, max_templates=max_templates)


def process_bam_parallel(qnamesortedbam, outh5, block_size=10000, threads=4, max_templates=None):

  p = Pool(threads)
  t_total = 0
  t0 = time.time()

  st = pd.HDFStore(outh5, mode='a', complevel=9, complib="blosc", format='t')

  t1 = time.time()
  for shard in p.imap_unordered(process_bam_section_w, get_bam_sections(qnamesortedbam, block_size, max_templates)):
    st.append('bdf', pd.DataFrame(shard), index=False)

    t1 = time.time()
    t_total += shard.shape[0]
    logger.debug('{} templates processed and written ({} templates/sec).'.format(t_total, t_total / (t1 - t0)))

  t2 = time.time()
  st.create_table_index('bdf', optlevel=9, kind='full')
  logger.debug('Took {} sec to index file'.format(t2 - t1))


def get_bam_sections(qnamesortedbam, block_size, max_templates=None):
  fp = pysam.AlignmentFile(qnamesortedbam)
  for cnt, _ in enumerate(fp):
    next(fp)
    next(fp)
    if cnt % block_size == 0:
      yield {
        'bam_name': qnamesortedbam,
        'file_offset': fp.tell(),
        'block_size': block_size
      }
    if max_templates is not None and cnt >= max_templates:
      break


def process_bam_section_w(args):
  """A thin wrapper to allow proper tracebacks when things go wrong in a thread

  :param args:
  :return:
  """
  import traceback
  try:
    return process_bam_section(args)
  except Exception as e:
    traceback.print_exc()
    print()
    raise e


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

def process_bam_section(args):

  df = np.empty((args['block_size'],), dtype=dfcols.get_bdf_cols().items())

  fp = pysam.AlignmentFile(args['bam_name'])
  fp.seek(args['file_offset'])

  for n in xrange(args['block_size']):
    r1 = next(fp, None)
    if r1 is None:
      df = df[:n]
      break

    r2 = next(fp)
    parse_read(r1, 'm1', df, n)
    parse_read(r2, 'm2', df, n)
    df['qname_hash'][n] = int(md5(r1.qname).hexdigest()[:8], 16)

  return df


def parse_read(read, mate, df, n):

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

  c_mapped = 0b10
  # Parse alignment
  long_insert = 0
  if read.is_unmapped:
    a_mapped = 0
    d = 1000000000
  else:
    a_mapped = 0b01
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

  # a_pos = ((read.reference_id + 1) << 29) | read.pos
  # c_pos = (chrom << 29) | pos

  df[mate + '_correct_pos'][n] = (chrom << 29) | pos
  df[mate + '_aligned_pos'][n] = ((read.reference_id + 1) << 29) | pos
  df[mate + '_mapped'][n] = a_mapped | c_mapped
  df[mate + '_d_error'][n] = d
  df[mate + '_MQ'][n] = read.mapping_quality

  # print(str(n) + ', ' + mate + ':' + str(d))

  tags = {k: v for k, v in read.get_tags()}
  for dtk in dfcols.graph_diagnostic_tags():
    df[mate + '_' + dtk[0]][n] = tags.get(dtk[0], 0)


if __name__ == '__main__':
  cli()