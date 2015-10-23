"""Functions to categorize reads for further analysis."""
import sqlite3 as sq
from itertools import izip
import array
import re

import h5py
import numpy as np
from pysam import AlignedSegment as pas

from mitty.lib.reads import old_style_cigar

import logging
logger = logging.getLogger(__name__)


def analyze_read(read, window=100, extended=False):
  """Given a read process the qname and read properties to determine the correct (CHROM, POS, CIGAR) and determine
  what kind of alignment errors were made on it

  :param read: a psyam AlignedSegment object
  :returns read_serial, chrom, cpy, ro, pos, cigar, chrom_c, pos_c, cigar_c, unmapped, t_start, t_end

  read_serial = read_serial * 10 + 0 or 1 (for mate1 or mate2 of read) for paired reads
  read_serial = read_serial for un-paired reads
  """
  early_exit_value = [None] * 15

  # Not counted
  if read.is_secondary:
    return early_exit_value

  ro_m, pos_m, rl_m, cigar_m = 0, 0, 0, ''  # These are the values passed in for unpaired reads
  # We should never actually fail this, unless a tool messes badly with the qname
  try:
    #  'read_serial|chrom|copy|ro|pos|rlen|cigar|ro|pos|rlen|cigar'
    if read.is_paired:
      if read.is_read1:
        rs, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m = read.qname.split('|')
      else:
        rs, chrom, cpy, ro_m, pos_m, rl_m, cigar_m, ro, pos, rl, cigar = read.qname.split('|')
      read_serial = int(rs) * 10 + (not read.is_read1)
    else:
      rs, chrom, cpy, ro, pos, rl, cigar = read.qname.split('|')[:9]
      read_serial = int(rs)
    ro, chrom, cpy, pos, rl, pos_m, rl_m = int(ro), int(chrom), int(cpy), int(pos), int(rl), int(pos_m), int(rl_m)
  except ValueError:
    logger.debug('Error processing qname: qname={:s}, chrom={:d}, pos={:d}'.format(read.qname, read.reference_id + 1, read.pos))
    return early_exit_value

  chrom_c, pos_c, cigar_c, unmapped = 1, 1, 1, 0

  if not extended:
    cigar = old_style_cigar(cigar)

  if read.is_unmapped:
    unmapped = 1
  else:
    if read.reference_id != chrom - 1:
      chrom_c, pos_c = 0, 0  # chrom wrong, so pos wrong too
    else:
      if check_read(read_pos=read.pos, read_cigar=read.cigarstring, correct_pos=pos, correct_cigar=cigar, window=window) != 0b000:
        pos_c = 0

    if read.cigarstring != cigar:  # TODO Use check read for this?
      cigar_c = 0

  return read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped


# TODO: Revise algorithm to properly work with reads inside long insertions
cigar_parser = re.compile(r'(\d+)(\D)')
def check_read(read_pos, read_cigar, correct_pos, correct_cigar, window):
  """

  :param read_pos:
  :param read_cigar:
  :param correct_pos:
  :param correct_cigar:
  :param window:
  :return: cat 3 bit value: category_code fragment indicating (pos bit, cigar bits)

  * For now the CIAGR bit is always set to 11 (fully correct) if pos is correct
  """
  cat = 0b111  # all wrong
  cigar_ops = cigar_parser.findall(correct_cigar)
  # Not comparing CIGARs right now. Will do for the future
  for cnt, op in cigar_ops:
    if op == '=' or op == 'M' or op == 'X':
      if -window <= read_pos - correct_pos <= window:
        cat = 0b000  # all correct
        break
      correct_pos += int(cnt)
    elif op == 'D':
      correct_pos += int(cnt)
  return cat


def count_reads_under_features(bam_fp, f_chrom_id, f_start, f_stop, f_chrom_cpy=None):
  """Count correct/total reads under given features

  :param bam_fp: perfect BAM file with read alignment analysis in extended tags
  :param f_chrom_id: chrom id of feature should allow us to fetch reads from the file
  :param f_chrom_cpy: copy of chromosome features come from. Set to none if this does not matter
  :param f_start: start of features -> Should be in ascending order
  :param f_stop: end of features
  :return: dict
   {
    "fully_outside_features": [correct, total],
    "templates_within_feature_but_read_outside": recarray same size as f_pos. Fields are correct, total
    "reads_within_feature": recarray same size as f_pos. Fields are correct, total
   }
  """
  non_feature_correct, non_feature_total = 0, 0
  templates_within_feature_but_read_outside = np.zeros(len(f_start), dtype=[('correct', 'uint32'), ('total', 'uint32')])
  twf_correct = templates_within_feature_but_read_outside['correct']
  twf_total = templates_within_feature_but_read_outside['total']
  reads_within_feature = np.zeros(len(f_start), dtype=[('correct', 'uint32'), ('total', 'uint32')])
  rwf_correct = reads_within_feature['correct']
  rwf_total = reads_within_feature['total']

  f_win_start = 0  # The window over the features
  f_win_stop = 0
  f_count = len(f_start)
  f_count_1 = f_count - 1

  for r in bam_fp.fetch(region=f_chrom_id):
    if f_chrom_cpy is not None:
      if r.get_tag('Zc') != f_chrom_cpy:  # Read comes from other chrom copy
        continue

    r_start = r.pos
    r_stop = r.get_tag('ZE')
    rm_start = r.pnext
    rm_stop = r.get_tag('Ze')

    template_start = min(r_start, rm_start)
    template_stop = max(r_stop, rm_stop)

    # Advance our variant window as needed
    if f_win_start < f_count_1:  # Only need to advance if there is something left
      while f_stop[f_win_start] < template_start:
        if f_win_start < f_count_1:
          f_win_start += 1
        else:
          break

    if f_win_stop < f_count_1:  # Only need to advance if there is something left
      while f_start[f_win_stop + 1] < template_stop:
        f_win_stop += 1
        if f_win_stop == f_count_1:
          break

    this_read_is_within_a_feature = False
    this_template_is_within_a_feature = False

    if f_win_start <= f_win_stop < f_count:
      # Now test this read against every feature
      for n in range(f_win_start, f_win_stop + 1):
        if r_start <= f_stop[n] and r_stop >= f_start[n]:
          this_read_is_within_a_feature = True
          rwf_total[n] += 1
          if r.get_tag('Xf') == 1: rwf_correct[n] += 1

      if not this_read_is_within_a_feature:
        for n in range(f_win_start, f_win_stop + 1):
          if rm_start <= f_stop[n] and rm_stop >= f_start[n]:
            this_template_is_within_a_feature = True
            twf_total[n] += 1
            if r.get_tag('Xf') == 1: twf_correct[n] += 1

    if not (this_read_is_within_a_feature or this_template_is_within_a_feature):
      non_feature_total += 1
      if r.get_tag('Xf') == 1: non_feature_correct += 1

  return {
    "fully_outside_features": [non_feature_correct, non_feature_total],
    "templates_within_feature_but_read_outside": templates_within_feature_but_read_outside,
    "reads_within_feature": reads_within_feature
  }


def find_nearest_variant(v1, v2):
  """For every variant in v1 find the nearest variant in v2 and note its distance and length. length = alt - ref such
  that SNP = 0, insertion > 0 and deletion < 0

  :param v1: list 1
  :param v2: list 2
  :return: numpy structured array the same length as v1 with fields 'dist' and 'length'
  """
  nearest_v2 = np.empty(v1.shape[0], dtype=[('dist', 'uint32'), ('length', 'int32')])
  nearest_dist, nearest_len = nearest_v2['dist'], nearest_v2['length']
  nearest_dist[:], nearest_len[:] = np.Inf, np.nan

  v2_size = v2.shape[0]
  v2_size_1 = v2_size - 1
  if v2_size == 0: return nearest_v2
  v2_idx = 0
  for n, v in enumerate(v1):
    d_min = abs(v2[v2_idx]['pos'] - v['pos'])
    while 1:
      if v2_idx < v2_size_1:
        d = abs(v2[v2_idx + 1]['pos'] - v['pos'])
        if d <= d_min:
          d_min = d
          v2_idx += 1
        else:
          break
      else:
        break
    nearest_dist[n] = d_min
    nearest_len[n] = len(v2[v2_idx]['alt']) - len(v2[v2_idx]['ref'])
  return nearest_v2


def categorize_read_counts_by_indel_length_and_nearest_variant(v1, v_read_counts, nearest_v2,
                                                               cat_counts=None,
                                                               max_v1_indel=100,
                                                               max_v2_indel=100,
                                                               max_dist=300):
  """v2 indels that are beyond the max range are piled into the first or last bins."""
  assert max_v1_indel > 0
  assert max_v2_indel > 0
  assert len(v1) == len(v_read_counts)
  assert len(v1) == len(nearest_v2)
  if cat_counts is None:
    cat_counts = {
      'counts': np.zeros((2 * max_v1_indel + 1, 2 * max_v2_indel + 1, max_dist + 1), dtype=[('correct', 'uint32'), ('total', 'uint32')]),
      'dimensions': [('indel_v1_size', range(-max_v1_indel, max_v1_indel + 1)),
                     ('indel_v2_size', range(-max_v1_indel, max_v1_indel + 1)),
                     ('dist', range(max_dist + 1))]
    }
  correct = cat_counts['counts']['correct']
  total = cat_counts['counts']['total']
  for n in xrange(v1.shape[0]):
    d_v1 = len(v1[n]['alt']) - len(v1[n]['ref'])
    #if abs(d_v1) > max_v1_indel: continue
    idx1 = max(0, min(2 * max_v1_indel, max_v1_indel + d_v1))
    idx2 = max(0, min(2 * max_v2_indel, nearest_v2[n]['length'] + max_v2_indel))
    idx3 = min(max_dist, nearest_v2[n]['dist'])
    total[idx1, idx2, idx3] += v_read_counts[n]['total']
    correct[idx1, idx2, idx3] += v_read_counts[n]['correct']
  return cat_counts