"""Associate reads in a BAM with variant calls.

Read through sorted BAM from aligner and place reads under FP and TP
If simulated reads, parse to fill out correct position

If available, read through inverted BAM, place reads under FN, fill out aligned pos

Finally, read through the qname sorted BAM and make sure all mates are present.

Remove any duplicates that are present


//samtools sort -n -@ 6 -o qname-run00001.S.pe.100x500.perfect.g-G1.t-0.8.17.bam unsorted-unindexed-run00001.S.pe.100x500.perfect.g-G1.t-0.8.17.bam

python -m cProfile -o myscript.cprof ~/Code/Mitty/mitty/benchmarking/reads2calls.py qname-run00001.S.pe.100x500.perfect.g-G1.t-0.8.17.bam eval.df.csv read-call.csv -v

"""
import logging
import time

import click
import pandas as pd
import pysam


logger = logging.getLogger(__name__)


@click.command()
@click.argument('qnamesortedbam')
@click.argument('evalvcfdf')
@click.argument('outcsv')
@click.option('--paired-reads/--single-end-reads', default=True, help='Are these paired end or single end reads')
@click.option('--simulated-reads/--real-reads', default=True, help='Are these simulated from Mitty? If so, parse qnames for truth positions')
@click.option('-v', count=True, help='Verbosity level')
def cli(qnamesortedbam, evalvcfdf, outcsv, paired_reads, simulated_reads, v):
  """Associate reads in a BAM with variant calls."""
  # print(alignerbam)
  # print(evalvcfdf)
  # print(outcsv)
  # print(paired_reads)
  # print(simulated_reads)

  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  logger.warning('This code only works for paired end simulated files ...')
  logger.warning('It is trvial to convert it to work for data with no qname information')

  #eval_df = pd.read_csv(evalvcfdf, dtype={'chrom': 'S10'}, compression='gzip' if evalvcfdf.endswith('gz') else None)
  #NOTE: Trying out only integers for chrom column
  eval_df = pd.read_csv(evalvcfdf, compression='gzip' if evalvcfdf.endswith('gz') else None)

  fp = pysam.AlignmentFile(qnamesortedbam)

  calls_over_bam_paired_simulated(
    eval_df, fp, progress_call_back=show_prog, progress_update_interval=100,
  ).to_csv(outcsv, index=False, compression='gzip' if outcsv.endswith('gz') else None)


def calls_over_bam_paired_simulated(eval_df, fp, progress_call_back=None, progress_update_interval=None):
  """For each read (pair) in the BAM, extract a_chrom, a_pos and c_chrom, c_pos if possible.
  Find the calls sitting above the aligned position and the correct position, if available and
  enter them into a dataframe. Return a list of dicts that can be put into a data frame.
  """
  seq_list = [s['SN'] for s in fp.header['SQ']]

  #from IPython import embed; embed()
  #exit()

  eval_df_mat = eval_df[['chrom', 'pos', 'pos_stop']].as_matrix()

  t0 = time.time()
  call_read_list = []
  for cnt, r1 in enumerate(fp):
    r2 = next(fp)
    new_calls = get_calls_over_template(eval_df, eval_df_mat, seq_list, r1, r2)
    if new_calls:
      call_read_list += new_calls

    if progress_call_back is not None:
      if cnt % progress_update_interval == 0:
        progress_call_back(cnt)
        if cnt > 500:
          break

  t1 = time.time()
  logger.debug('{} templates in {} sec ({} templates/sec)'.format(cnt, t1 - t0, cnt/(t1 - t0)))

  # All this to organize the columns how we like them
  df = pd.DataFrame(call_read_list)
  original_columns = df.columns.tolist()
  strictly_ordered_columns = eval_df.columns.tolist() + [
    'under_feature', 'd_error', 'qname', 'mate',
    'c_chrom', 'c_pos', 'c_cigar', 'long_insert', 'a_chrom', 'a_pos', 'a_cigar',
    'MQ', 'mapped'
  ]
  for k in strictly_ordered_columns:
    original_columns.remove(k)
  new_columns = strictly_ordered_columns + original_columns

  return pd.DataFrame(df, columns=new_columns).sort_values(columns=['chrom', 'pos'])


def get_calls_over_template(eval_df, eval_df_mat, seq_list, r1, r2):
  """given a read, find the calls that go over it"""

  call_read_list = []

  r1_d = parse_qname(r1, seq_list)
  r2_d = parse_qname(r2, seq_list)

  r1_a_ch, r1_a_p0, r1_a_p1 = r1_d['a_chrom'], r1_d['a_pos'], r1_d['a_pos'] + r1.rlen  # Location aligned
  r1_c_ch, r1_c_p0, r1_c_p1 = r1_d['c_chrom'], r1_d['c_pos'], r1_d['c_pos'] + r1.rlen  # Location should have been aligned

  r2_a_ch, r2_a_p0, r2_a_p1 = r2_d['a_chrom'], r2_d['a_pos'], r2_d['a_pos'] + r2.rlen  # Location aligned
  r2_c_ch, r2_c_p0, r2_c_p1 = r2_d['c_chrom'], r2_d['c_pos'], r2_d['c_pos'] + r2.rlen  # Location should have been aligned

  # All calls that are over where the read is aligned, or should have been aligned
  # r1_aligned_over = (eval_df['chrom'] == r1_a_ch) & (r1_a_p0 < eval_df['pos']) & (r1_a_p1 > eval_df['pos_stop'])
  # r1_correct_over = (eval_df['chrom'] == r1_c_ch) & (r1_c_p0 < eval_df['pos']) & (r1_c_p1 > eval_df['pos_stop'])
  #
  # r2_aligned_over = (eval_df['chrom'] == r2_a_ch) & (r2_a_p0 < eval_df['pos']) & (r2_a_p1 > eval_df['pos_stop'])
  # r2_correct_over = (eval_df['chrom'] == r2_c_ch) & (r2_c_p0 < eval_df['pos']) & (r2_c_p1 > eval_df['pos_stop'])

  # r1_aligned_over = (eval_df_ch == r1_a_ch) & (r1_a_p0 < eval_df_pos) & (r1_a_p1 > eval_df_pos_stop)
  # r1_correct_over = (eval_df_ch == r1_c_ch) & (r1_c_p0 < eval_df_pos) & (r1_c_p1 > eval_df_pos_stop)
  #
  # r2_aligned_over = (eval_df_ch == r2_a_ch) & (r2_a_p0 < eval_df_pos) & (r2_a_p1 > eval_df_pos_stop)
  # r2_correct_over = (eval_df_ch == r2_c_ch) & (r2_c_p0 < eval_df_pos) & (r2_c_p1 > eval_df_pos_stop)

  r1_aligned_over = (eval_df_mat[:,0] == r1_a_ch) & (r1_a_p0 < eval_df_mat[:,1]) & (r1_a_p1 > eval_df_mat[:,2])
  r1_correct_over = (eval_df_mat[:,0] == r1_c_ch) & (r1_c_p0 < eval_df_mat[:,1]) & (r1_c_p1 > eval_df_mat[:,2])
  
  r2_aligned_over = (eval_df_mat[:,0] == r2_a_ch) & (r2_a_p0 < eval_df_mat[:,1]) & (r2_a_p1 > eval_df_mat[:,2])
  r2_correct_over = (eval_df_mat[:,0] == r2_c_ch) & (r2_c_p0 < eval_df_mat[:,1]) & (r2_c_p1 > eval_df_mat[:,2])

  for call in eval_df[(r1_aligned_over | r1_correct_over) & ~(r2_aligned_over | r2_correct_over)].iterrows():
    c_dict = call[1].to_dict()
    call_read_list.append(
      {k: v for k, v in c_dict.items() + r1_d.items() + [('under_feature', 1)]}
    )
    call_read_list.append(
      {k: v for k, v in c_dict.items() + r2_d.items() + [('under_feature', 0)]}
    )

  for call in eval_df[~(r1_aligned_over | r1_correct_over) & (r2_aligned_over | r2_correct_over)].iterrows():
    c_dict = call[1].to_dict()
    call_read_list.append(
      {k: v for k, v in c_dict.items() + r1_d.items() + [('under_feature', 0)]}
    )
    call_read_list.append(
      {k: v for k, v in c_dict.items() + r2_d.items() + [('under_feature', 1)]}
    )

  for call in eval_df[(r1_aligned_over | r1_correct_over) & (r2_aligned_over | r2_correct_over)].iterrows():
    c_dict = call[1].to_dict()
    call_read_list.append(
      {k: v for k, v in c_dict.items() + r1_d.items() + [('under_feature', 1)]}
    )
    call_read_list.append(
      {k: v for k, v in c_dict.items() + r2_d.items() + [('under_feature', 1)]}
    )

  return call_read_list


def parse_qname(read, seq_list):

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

  # Parse alignment
  long_insert = 0
  if read.is_unmapped:
    mapped = 0
    d = 1000000000
  else:
    mapped = 1
    if read.reference_id != chrom - 1:
      d = 1000000000
    else:  # Analyze the correctness by checking each breakpoint
      # Corner case, our special cigar for indicating reads inside an insertion
      # We use S or I for this
      if cigar_ops[0][0] in [1, 4] and len(cigar_ops) == 1:  # S,I
        d = read.pos - pos
        long_insert = 1
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

    # This is a very strict checking of the CIGAR
    # TODO: Make this more clever
    if read.cigarstring != cigar:
      cigar_c = 0

  r_dict = {
    'd_error': d,  # Fill this out with proper metric
    'qname': read.qname,
    'mate': 0 if read.is_read1 else 1,
    'a_chrom': read.reference_id + 1, #NOTE: trying out ints for chrom# read.reference_name,
    'a_pos': read.pos,
    'a_cigar': a_cigar,
    'c_chrom': chrom,  #NOTE: trying out ints for chrom#  seq_list[chrom - 1],
    'c_pos': int(pos),
    'c_cigar': cigar,
    'long_insert': long_insert,
    'MQ': read.mapping_quality,
    'mapped': mapped
  }
  for k, v in read.get_tags():
    r_dict[k] = v
  return r_dict


def show_prog(cnt):
  logger.debug('{} read pairs'.format(cnt))


if __name__ == '__main__':
  cli()