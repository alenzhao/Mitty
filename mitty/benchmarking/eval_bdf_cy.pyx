from cpython cimport array
import array
import time
import logging

cimport cython
import numpy as np
cimport numpy as np


logger = logging.getLogger(__name__)

# Pure python version Took 19.30675101s to process 60000/244546 reads/calls (3107.72123 reads/s) with 2 threads
# Cython with no types Took 14.28495002s to process 60000/244546 reads/calls (4200.224707 reads/s) with 2 threads

def get_templates_over_calls(bam_df, eval_df, mate='m1', mode='a'):
  """

  :param bam_df: Whole or part of a bdf
  :param eval_df: Whole (or part) of an evdf
  :param mate: m1 or m2 : sort on mate1 or mate2
  :param mode: 'a' or 'c' : sort on aligned or correct position
  :return:
  """
  t0 = time.time()
  call_read_rows = []

  feature_key = [('{}_{}_under_feature'.format(mate, mode), 1)]

  p1_key = '{}_{}_p1'.format(mate, mode)  # e.g. 'm1_a_p1'
  p2_key = '{}_{}_p2'.format(mate, mode)  # e.g. 'm1_a_p2'

  d_sorted = bam_df.sort_values(by=[p1_key])

  r_p1 = d_sorted[p1_key].values
  r_p2 = d_sorted[p2_key].values

  v_p1 = eval_df['call_p1'].values
  v_p2 = eval_df['call_p2'].values

  call_cnt = len(eval_df)
  t_idx1, t_idx2, template_cnt = 0, 0, len(bam_df)

  for call_index in range(call_cnt):
    if t_idx1 >= template_cnt:
      break

    # Move read window such that all reads in window are over call
    while t_idx1 < template_cnt and (r_p2[t_idx1] < v_p1[call_index]):  # Move upper edge of window
      t_idx1 += 1

    while t_idx2 < template_cnt and (r_p1[t_idx2] < v_p2[call_index]):  # Move lower edge of window
      t_idx2 += 1
      # There are some shenanigans here with the last read in the file
      # (we don't consider it) but we don't care.

    if t_idx2 > t_idx1:
      citems = eval_df.iloc[call_index].to_dict().items()
      for r in d_sorted.iloc[t_idx1:t_idx2].iterrows():
        call_read_rows.append(
          {k: v for k, v in citems + r[1].to_dict().items() + feature_key}
        )

  t1 = time.time()
  logger.debug('Matched {} reads against {} calls in {:0.3} s'.format(template_cnt, call_cnt, t1 - t0))

  return call_read_rows
