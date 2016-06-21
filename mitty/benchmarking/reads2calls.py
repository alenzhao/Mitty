"""Associate reads in a BAM with variant calls.

Read through sorted BAM from aligner and place reads under FP and TP
If simulated reads, parse to fill out correct position

If available, read through inverted BAM, place reads under FN, fill out aligned pos

Finally, read through the qname sorted BAM and make sure all mates are present.

Remove any duplicates that are present



python ~/Code/Mitty/mitty/benchmarking/reads2calls.py unsorted-unindexed-run00001.S.pe.100x500.perfect.g-G1.t-0.8.17.bam eval.df.csv read-call.vcf -p

"""
import logging
import io
import os

import click
import pandas as pd
import pysam


logger = logging.getLogger(__name__)


@click.command()
@click.argument('alignerbam')
@click.argument('evalvcfdf')
@click.argument('outcsv')
@click.option('--paired-reads/--single-end-reads', default=True, help='Are these paired end or single end reads')
@click.option('--simulated-reads/--real-reads', default=True, help='Are these simulated from Mitty? If so, parse qnames for truth positions')
@click.option('-p', is_flag=True, help='Show progressbar')
@click.option('-v', count=True, help='Verbosity level')
def cli(alignerbam, evalvcfdf, outcsv, paired_reads, simulated_reads, p, v):
  """Associate reads in a BAM with variant calls."""
  # print(alignerbam)
  # print(evalvcfdf)
  # print(outcsv)
  # print(paired_reads)
  # print(simulated_reads)

  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  eval_df = pd.read_csv(evalvcfdf, dtype={'chrom': 'S10'}, compression='gzip' if evalvcfdf.endswith('gz') else None)

  fp = pysam.AlignmentFile(alignerbam)
  total_file_size = os.path.getsize(alignerbam)
  progress_bar_update_interval = 0.01 * total_file_size

  with click.progressbar(
    length=total_file_size, label='Processing BAM', file=None if p else io.BytesIO()) as bar:
    calls_over_bam_paired_simulated(
      eval_df, fp, progress_call_back=bar, progress_bar_update_interval=progress_bar_update_interval
    ).to_csv(outcsv, index=False, compression='gzip' if outcsv.endswith('gz') else None)


def calls_over_bam_paired_simulated(eval_df, fp, progress_call_back=None, progress_bar_update_interval=None):
  """For each read (pair) in the BAM, extract a_chrom, a_pos and c_chrom, c_pos if possible.
  Find the calls sitting above the aligned position and the correct position, if available and
  enter them into a dataframe. Return a list of dicts that can be put into a data frame.
  """
  logger.warning('This code only works for paired end simulated files ...')
  logger.warning('It is trvial to convert it to work for data with no qname information')
  logger.warning('It is hardcoded for speed purposes')

  # fp = pysam.AlignmentFile(alignerbam)
  # progress_bar_update_interval = 0.01 * (fp.mapped + fp.unmapped) / 2
  seq_list = [s['SN'] for s in fp.header['SQ']]

  call_read_list = []
  last_file_pos = 0
  for r1 in fp:
    r2 = next(fp)
    new_calls = get_calls_over_template(eval_df, seq_list, r1, r2)
    if new_calls:
      call_read_list += new_calls

    if progress_call_back is not None:
      if fp.tell() - last_file_pos > progress_bar_update_interval:
        progress_call_back.update(fp.tell() - last_file_pos)
        last_file_pos = fp.tell()
      if fp.tell() > 4 * progress_bar_update_interval:
        break

  #  from IPython import embed; embed()

  # All this to organize the columns how we like them
  df = pd.DataFrame(call_read_list)
  original_columns = df.columns.tolist()
  strictly_ordered_columns = eval_df.columns.tolist() + [
    'under_feature', 'd_error', 'qname', 'mate',
    'c_chrom', 'c_pos', 'c_cigar', 'a_chrom', 'a_pos', 'a_cigar',
    'MQ', 'mapped'
  ]
  for k in strictly_ordered_columns:
    original_columns.remove(k)
  new_columns = strictly_ordered_columns + original_columns

  return pd.DataFrame(df, columns=new_columns)


def get_calls_over_template(eval_df, seq_list, r1, r2):
  """given a read, find the calls that go over it"""

  call_read_list = []

  r1_d = parse_qname(r1, seq_list)
  r2_d = parse_qname(r2, seq_list)

  r1_a_ch, r1_a_p0, r1_a_p1 = r1_d['a_chrom'], r1_d['a_pos'], r1_d['a_pos'] + r1.rlen  # Location aligned
  r1_c_ch, r1_c_p0, r1_c_p1 = r1_d['c_chrom'], r1_d['c_pos'], r1_d['c_pos'] + r1.rlen  # Location should have been aligned

  r2_a_ch, r2_a_p0, r2_a_p1 = r2_d['a_chrom'], r2_d['a_pos'], r2_d['a_pos'] + r2.rlen  # Location aligned
  r2_c_ch, r2_c_p0, r2_c_p1 = r2_d['c_chrom'], r2_d['c_pos'], r2_d['c_pos'] + r2.rlen  # Location should have been aligned

  # All calls that are over where the read is aligned, or should have been aligned
  r1_aligned_over = (eval_df['chrom'] == r1_a_ch) & (r1_a_p0 < eval_df['pos']) & (r1_a_p1 > eval_df['pos_stop'])
  r1_correct_over = (eval_df['chrom'] == r1_c_ch) & (r1_c_p0 < eval_df['pos']) & (r1_c_p1 > eval_df['pos_stop'])

  r2_aligned_over = (eval_df['chrom'] == r2_a_ch) & (r2_a_p0 < eval_df['pos']) & (r2_a_p1 > eval_df['pos_stop'])
  r2_correct_over = (eval_df['chrom'] == r2_c_ch) & (r2_c_p0 < eval_df['pos']) & (r2_c_p1 > eval_df['pos_stop'])

  for call in eval_df.loc[(r1_aligned_over | r1_correct_over) & ~(r2_aligned_over | r2_correct_over)].iterrows():
    c_dict = call[1].to_dict()
    call_read_list.append(
      {k: v for k, v in c_dict.items() + r1_d.items() + [('under_feature', 1)]}
    )
    call_read_list.append(
      {k: v for k, v in c_dict.items() + r2_d.items() + [('under_feature', 0)]}
    )

  for call in eval_df.loc[~(r1_aligned_over | r1_correct_over) & (r2_aligned_over | r2_correct_over)].iterrows():
    c_dict = call[1].to_dict()
    call_read_list.append(
      {k: v for k, v in c_dict.items() + r1_d.items() + [('under_feature', 0)]}
    )
    call_read_list.append(
      {k: v for k, v in c_dict.items() + r2_d.items() + [('under_feature', 1)]}
    )

  for call in eval_df.loc[(r1_aligned_over | r1_correct_over) & (r2_aligned_over | r2_correct_over)].iterrows():
    c_dict = call[1].to_dict()
    call_read_list.append(
      {k: v for k, v in c_dict.items() + r1_d.items() + [('under_feature', 1)]}
    )
    call_read_list.append(
      {k: v for k, v in c_dict.items() + r2_d.items() + [('under_feature', 1)]}
    )

  return call_read_list


def parse_qname(read, seq_list):
  if read.is_paired:
    if read.is_read1:
      rs, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m = read.qname.split('|')
    else:
      rs, chrom, cpy, ro_m, pos_m, rl_m, cigar_m, ro, pos, rl, cigar = read.qname.split('|')
  else:
    rs, chrom, cpy, ro, pos, rl, cigar = read.qname.split('|')[:7]

  r_dict = {
    'd_error': 0,  # Fill this out with proper metric
    'qname': read.qname,
    'mate': 0 if read.is_read1 else 1,
    'a_chrom': read.reference_name,
    'a_pos': read.pos,
    'a_cigar': read.cigarstring,
    'c_chrom': seq_list[int(chrom) - 1],
    'c_pos': int(pos),
    'c_cigar': cigar,
    'MQ': read.mapping_quality,
    'mapped': 1 - read.is_unmapped
  }
  for k, v in read.get_tags():
    r_dict[k] = v
  return r_dict


if __name__ == '__main__':
  cli()