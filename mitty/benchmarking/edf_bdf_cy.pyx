import time
import logging

cimport cython
import numpy as np
cimport numpy as np

import pandas as pd


logger = logging.getLogger(__name__)



# Pure python version Took 19.30675101s to process 60000/244546 reads/calls (3107.72123 reads/s) with 2 threads
# Cython with no types Took 14.28495002s to process 60000/244546 reads/calls (4200.224707 reads/s) with 2 threads
# Cython with types and better indexing  Took 5024.602472s to process 30000000/244546 templates/calls (5970.62159 templates/s) with 2 threads
def get_templates_over_calls(bam_df, eval_df):
  """

  :param bam_df: Whole or part of a bdf
  :param eval_df: Whole (or part) of an evdf
  :return:
  """
  t0 = time.time()
  call_read_rows = []

  cdef:
    np.ndarray[long, ndim=1] c_p1, c_p2
    np.ndarray[long, ndim=2] r_p1, r_p2
    np.ndarray[size_t, ndim=2] r_srt_idx

    size_t call_index = 0, call_cnt = len(eval_df), template_cnt = len(bam_df), idx
    np.ndarray[size_t, ndim=1] t_idx1 = np.zeros(4, dtype=np.uintp), t_idx2 = np.zeros(4, dtype=np.uintp)
    int condition = 0

  # c_p1, c_p2 are the start and stop positions of the calls

  c_p1 = eval_df['call_p1'].values
  c_p2 = eval_df['call_p2'].values

  # r_p1, r_p2 are the (unsorted) arrays of read start and stop positions for combinations of
  # mate and aligned/correct

  # r_srt_idx are the indexes that sort each axis of r_p1 in ascending order
  # These are what we use to find the reads under each call

  r_p1 = np.empty((bam_df.shape[0], 4), dtype=long)
  r_p2 = np.empty((bam_df.shape[0], 4), dtype=long)
  r_srt_idx = np.empty((bam_df.shape[0], 4), dtype=np.uintp)


  flag_headings = {}
  flag_meanings = []
  cntr = 0
  for mate in ['m1', 'm2']:
    for mode in ['a', 'c']:
      p1_key = '{}_{}_p1'.format(mate, mode)  # e.g. 'm1_a_p1'
      p2_key = '{}_{}_p2'.format(mate, mode)  # e.g. 'm1_a_p2'
      r_p1[:, cntr] = bam_df[p1_key].values
      r_p2[:, cntr] = bam_df[p2_key].values
      r_srt_idx[:, cntr] = np.argsort(r_p1[:, cntr])
      flag_meanings.append('{}_{}_here'.format(mate, 'aligned' if mode == 'a' else 'should_be'))
      flag_headings['{}_{}_here'.format(mate, 'aligned' if mode == 'a' else 'should_be')] = 0
      cntr += 1
  read_type = pd.Series(flag_headings, dtype=int)


  for call_index in range(call_cnt):
    advance_indexes(
      call_index, c_p1, c_p2,
      r_p1, r_p2, r_srt_idx,
      t_idx1, t_idx2, template_cnt)
    call_read_rows += prepare_rows(
      eval_df, bam_df, call_index, r_srt_idx, t_idx1, t_idx2, read_type, flag_meanings)

  t1 = time.time()
  logger.debug('Matched {} templates against {} calls in {:0.3} s'.format(template_cnt, call_cnt, t1 - t0))

  return call_read_rows


cdef advance_indexes(
  size_t call_index, np.ndarray[long, ndim=1] c_p1, np.ndarray[long, ndim=1] c_p2,
  np.ndarray[long, ndim=2] r_p1, np.ndarray[long, ndim=2] r_p2,
  np.ndarray[size_t, ndim=2] r_srt_idx,
  np.ndarray[size_t, ndim=1] t_idx1, np.ndarray[size_t, ndim=1] t_idx2, size_t template_cnt):

  for condition in range(4):
    advance_index(
      call_index, c_p1, c_p2,
      r_p1, r_p2, r_srt_idx,
      t_idx1, t_idx2, template_cnt, condition)

cdef advance_index(size_t call_index, np.ndarray[long, ndim=1] c_p1, np.ndarray[long, ndim=1] c_p2,
                   np.ndarray[long, ndim=2] r_p1, np.ndarray[long, ndim=2] r_p2,
                   np.ndarray[size_t, ndim=2] r_srt_idx,
                   np.ndarray[size_t, ndim=1] t_idx1, np.ndarray[size_t, ndim=1] t_idx2, size_t template_cnt,
                   int condition):
  if t_idx1[condition] >= template_cnt:
    return

  # Move read window such that all reads in window are over call
  while t_idx1[condition] < template_cnt and (r_p2[r_srt_idx[t_idx1[condition], condition], condition] < c_p1[call_index]):  # Move upper edge of window
    t_idx1[condition] += 1

  while t_idx2[condition] < template_cnt and (r_p1[r_srt_idx[t_idx2[condition], condition], condition] < c_p2[call_index]):  # Move lower edge of window
    t_idx2[condition] += 1
    # There are some shenanigans here with the last read in the file
    # (we don't consider it) but we don't care.


cdef prepare_rows(
  eval_df, bam_df, size_t call_index,
  np.ndarray[size_t, ndim=2] r_srt_idx,
  np.ndarray[size_t, ndim=1] t_idx1, np.ndarray[size_t, ndim=1] t_idx2, read_type, flag_meanings):
  cdef:
    size_t condition, idx, last_qname_idx
    np.ndarray[size_t, ndim=2] qname_idx_flag  # read index and flag

  rows = []

  npa = [
    [idx, condition]
    for condition in range(4)
    for idx in r_srt_idx[t_idx1[condition]:t_idx2[condition], condition]
  ]

  if len(npa) > 0:
    qname_idx_flag = np.array(npa, dtype=np.uintp)
    qname_idx_flag = qname_idx_flag[np.argsort(qname_idx_flag[:, 0]), :]

    call_s = eval_df.ix[call_index]
    last_qname_idx = qname_idx_flag[0, 0]
    this_read_type = pd.Series(read_type)
    for idx in range(qname_idx_flag.shape[0]):
      if qname_idx_flag[idx, 0] == last_qname_idx:
        this_read_type[flag_meanings[qname_idx_flag[idx, 1]]] = 1
      else:
        rows += [pd.concat((call_s, bam_df.iloc[last_qname_idx], this_read_type))]

        last_qname_idx = <int>qname_idx_flag[idx, 0]
        this_read_type = pd.Series(read_type)
        this_read_type[flag_meanings[qname_idx_flag[idx, 1]]] = 1

    rows += [pd.concat((call_s, bam_df.iloc[last_qname_idx], this_read_type))]

  return rows