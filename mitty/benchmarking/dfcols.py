"""
A lot of the data prep ivolves figuring out what columns to save, what datatype to use, what to index.
It's all here in one place with a bit of commentary"""
from hashlib import md5
from collections import OrderedDict as od

import numpy as np


def call_hash(row):
  return int(md5(','.join([row['#CHROM'], str(row['POS']), row['REF'], row['ALT']])).hexdigest()[:8], 16)


def call_type_encoding(row):
  if row['truth'] == 'TP':
    return 0
  elif row['truth'] == 'FN':
    return 1
  else:
    return 2  # FP


def get_edf_data_cols():
  """Note that the evcf is a CSV and does not really have data columns defined."""
  return od([
    # These are computed columns that are useful for the CR structure
    ('call_hash', np.uint64),
    ('call_type', np.uint8),     # Computed from truth/query columns
    ('call_p1', np.uint64),      # (chrom << 29) | pos
    ('call_p2', np.uint64),      # (chrom << 29) | (pos + variant_size)
    ('variant_size', np.int32),  # len(alt) - len(ref)
  ])


def get_edf_cols(chrom_id_fmt=None):
  """
  EDF is basically a copy of the eval.vcf with columns cleaned up and parsed to be more easily processed.
  This is a gzipped csv file, so we don't really need data types. A subset of these columns go into the
  cr (calls/reads) table which is HDF5 and it needs types then.
  """
  dc = get_edf_data_cols()  # These are computed columns that are useful for the CR structure

  return od(dc.items() + [
      # These are the original VCF columns
    ('call_chrom', chrom_id_fmt or 'a10'),  # The evcf converter may modify this to fit the longest chrom id string
    ('call_pos', np.uint32),
    ('ref', object),
    ('alt', object),
    ('ROC_thresh', np.float16),
    ('truth', 'a2'),  # FP, FN, TP
    ('truthGT', 'a5'),
    ('query', 'a2'),  # FP, FN, TP
    ('queryGT', 'a5')
    ])


def get_bdf_cols():
  """
  BDF is a subset of BAM read information, suitable for aligner analysis summaries and for computing
  the CR structure. Not suitable for recreating sample FASTQ files - we would use the qname_hashes and
  stream through the original BAM file for that.

  The following particular numerical information for each read is processed and stored in an HDF file.
  (We choose an HDF5 rather than a gzipped csv so that we can efficiently partition the data in the file
  and process it in parallel. A CSV file would require a central reader that passes read chunks to
  child processes)

  We index all the columns.

  """
  qh = [('qname_hash', np.uint64)]  # allows us to find reads by qname. Immune to sort order etc.
  # Regular info for simulated BAMs. Needs to be duplicated for mate1 and mate2
  reg = [
    ('correct_pos', np.uint64),  # these are computed as (chrom << 29) | pos
    ('aligned_pos', np.uint64),
    ('mapped', np.uint8),        # 2bits,   b0 = 0,1 aligned is mapped or not b1 = 0,1 correct is mapped or not
    ('d_error', np.int64),
    ('MQ', np.uint8)
  ]
  graph = graph_diagnostic_tags()  # Extended tags from the graph aligner

  return od(qh + [(m + '_' + tag[0], tag[1]) for tag in reg + graph for m in ['m1', 'm2']])


def graph_diagnostic_tags():
  return [
    # The following are extended tags from the graph aligner
    ('YS', np.uint8),
    ('YA', np.uint8),
    ('Ym', np.uint16),
    ('UQ', np.uint8),
    ('YI', np.uint16),
    ('YQ', np.uint16)
  ]


#
# def get_cr_cols(data_only=False):
#   """
#   CR - associates reads with calls
#   :return:
#   """
#   if data_only:
#     return get_edf_data_cols()
#
#   return od(get_edf_data_cols().items() + get_bdf_cols().items())
#
#
# read_info_cols = [
#   'd_error',   # Alignment error metric
#   'a_mapped',  # Aligned mapped status
#   'a_p1',      # Aligned pos = chrom << 29 | pos
#   'a_p2',      # End pos of read
#   'a_cigar',
#   'c_mapped',  # Correct mapped status
#   'c_p1',      # Correct pos
#   'c_p2',      # End pos of read
#   'c_cigar',
#   'MQ'
# ]
#
# gral_tag_cols = [
# #  'XG',
#   'YG',
#   'Yg',
#   'XI',
# #  'XB',
#   'XE',
#   'YX',
#   'Yx',
#   'YS',
#   'YA',
#   'YQ',
#   'Ym',
#   'UQ'
# ]
#
# # These are the columns the bam2df file will have
# bdf_cols = ['qname', 'qname_hash'] + [m + t for m in ['m1_', 'm2_'] for t in read_info_cols + gral_tag_cols]
#
#
# # For the edf-bdf structure
# # These are the cols from the call columns we will take
# edf_bdf_call_cols = [
#   'call_hash',
#   'ROC_thresh',  # Is this worth it?
#   'call_type',
#   'variant_size'
# ]
# # And these are the cols from the read we will take
# r_cols = ['d_error', 'MQ', 'UQ', 'YQ', 'YS', 'YA', 'Ym']
# edf_bdf_read_cols = ['qname_hash'] + [m + t for m in ['m1_', 'm2_'] for t in r_cols] + ['align_flag']
#
#
# # For now, not using a multiple table store, for simplicity, may revisit if things get slow/clunky
# # Just trying out restricting the indexes
# edf_bdf_data_cols = [
#   'call_hash',
#   'call_type',
#   'variant_size'
# ]
