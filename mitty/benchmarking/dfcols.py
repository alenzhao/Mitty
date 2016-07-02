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
    ('correct_p1', np.uint64),  # these are computed as (chrom << 29) | pos
    ('correct_p2', np.uint64),
    ('aligned_p1', np.uint64),
    ('aligned_p2', np.uint64),
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