"""Just a set of lists indicating our preferred order for columns for some standard data frames."""

evdf_cols = [
  u'call_chrom',
  u'call_pos',
  u'ref',
  u'alt',
  u'ROC_thresh',
  u'truth',
  u'truthGT',
  u'query',
  u'queryGT',
  u'variant_size',
  u'call_p1',
  u'call_p2'
]


read_info_cols = [
  'd_error',   # Alignment error metric
  'a_mapped',  # Aligned mapped status
  'a_p1',      # Aligned pos = chrom << 29 | pos
  'a_p2',      # End pos of read
  'a_cigar',
  'c_mapped',  # Correct mapped status
  'c_p1',      # Correct pos
  'c_p2',      # End pos of read
  'c_cigar',
  'MQ'
]

gral_tag_cols = [
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

bdf_cols = ['qname'] + [m + t for m in ['m1_', 'm2_'] for t in read_info_cols + gral_tag_cols]

calls_reads_cols = evdf_cols + bdf_cols + ['{}_{}_under_feature'.format(mate, mode) for mate in ['m1', 'm2'] for mode in ['a', 'c']]
