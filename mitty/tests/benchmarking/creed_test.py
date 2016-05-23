from nose.tools import assert_raises

import pysam

import mitty.benchmarking.creed as creed


# class MyRead:
#   def __init__(self, qname, secondary, paired, read1, unmapped, reference_id, pos, cigarstring):
#     self.qname, self.is_secondary, self.is_paired, self.is_read1 = qname, secondary, paired, read1
#     self.is_unmapped, self.reference_id, self.pos, self.cigarstring = unmapped, reference_id, pos, cigarstring


def analyze_read_test0():
  """analyze read: ignore secondary alignment"""
  r = pysam.AlignedSegment()
  r.is_secondary = True
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d, long_insert, long_insert_m = creed.analyze_read(r)

  assert read_serial is None


def analyze_read_test1():
  """analyze read: ignore supplementary alignment"""
  r = pysam.AlignedSegment()
  r.is_supplementary = True
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d, long_insert, long_insert_m = creed.analyze_read(r)

  assert read_serial is None


def analyze_read_test2():
  """analyze read: raise runtime error for malformed qname"""
  r = pysam.AlignedSegment()
  r.qname = 'To forgive is divine'

  assert_raises(RuntimeError, creed.analyze_read, r)


# @read_serial|chrom|copy|ro|pos|rlen|cigar|ro|pos|rlen|cigar

def analyze_read_test3():
  """analyze read: parse SE qname"""
  r = pysam.AlignedSegment()
  r.qname = '1|1|0|0|1|100|100='
  r.is_paired = False
  r.cigarstring = '100='
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d, long_insert, long_insert_m = creed.analyze_read(r, extended=True)

  assert read_serial == 1
  assert chrom == 1
  assert cpy == 0
  assert ro == 0
  assert pos == 1
  assert rl == 100
  assert cigar == '100='


def analyze_read_test4():
  """analyze read: extended/traditional CIGAR conversion"""
  r = pysam.AlignedSegment()
  r.qname = '1|1|0|0|1|100|100='
  r.is_paired = False
  r.cigarstring = '100='
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d, long_insert, long_insert_m = creed.analyze_read(r)

  assert cigar == '100M'


def analyze_read_test5():
  """analyze read: parse paired qname"""
  r = pysam.AlignedSegment()
  r.qname = '1|1|0|0|1|100|100=|1|201|100|100='
  r.is_paired = True
  r.is_read1 = True
  r.cigarstring = '100='
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d, long_insert, long_insert_m = creed.analyze_read(r)

  assert read_serial == 10
  assert chrom == 1
  assert cpy == 0
  assert ro == 0
  assert pos == 1
  assert rl == 100
  assert cigar == '100M'

  assert ro_m == 1
  assert pos_m == 201
  assert rl_m == 100
  assert cigar_m == '100M'


def analyze_read_test6():
  """analyze read: check reference read"""
  r = pysam.AlignedSegment()
  r.qname = '1|1|0|0|1|100|100=|1|201|100|100='
  r.is_paired = True

  # Read 1
  r.is_read1 = True
  r.cigarstring = '100M'
  r.reference_id = 0
  r.pos = 1
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d, long_insert, long_insert_m = creed.analyze_read(r)

  assert chrom_c == 1
  assert pos_c == 1
  assert cigar_c == 1
  assert d == 0

  # Read 2
  r.is_read1 = False
  r.cigarstring = '100M'
  r.reference_id = 0
  r.pos = 201
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d, long_insert, long_insert_m = creed.analyze_read(r)

  assert chrom_c == 1
  assert pos_c == 1
  assert cigar_c == 1
  assert d == 0

  # Read 1 off by 10
  r.is_read1 = True
  r.cigarstring = '1M10D99M'
  r.reference_id = 0
  r.pos = 11
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d, long_insert, long_insert_m = creed.analyze_read(r)

  assert chrom_c == 1
  assert pos_c == 0
  assert cigar_c == 0
  assert d == 10

  # Read 1 is unmapped
  r.is_read1 = True
  r.is_unmapped = True
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d, long_insert, long_insert_m = creed.analyze_read(r)

  assert chrom_c == 0
  assert pos_c == 0
  assert cigar_c == 0
  assert d == creed.MAX_D_ERROR + 1

  # Read 1 on wrong chrom
  r.is_unmapped = False
  r.is_read1 = True
  r.cigarstring = '100M'
  r.reference_id = 1
  r.pos = 1
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d, long_insert, long_insert_m = creed.analyze_read(r)

  assert chrom_c == 0
  assert pos_c == 0
  assert cigar_c == 1
  assert d == creed.MAX_D_ERROR + 1


def analyze_read_test7():
  """analyze read: check read with insertion"""
  r = pysam.AlignedSegment()
  r.qname = '1|1|0|0|21|100|20S80=|1|201|100|100='
  r.is_paired = True
  r.is_read1 = True

  r.cigarstring = '100M'
  r.reference_id = 0
  r.pos = 21
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d, long_insert, long_insert_m = creed.analyze_read(r)

  assert chrom_c == 1
  assert pos_c == 1
  assert cigar_c == 0
  assert d == 0

  r.cigarstring = '20S80M'
  r.reference_id = 0
  r.pos = 10
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d, long_insert, long_insert_m = creed.analyze_read(r)

  assert chrom_c == 1
  assert pos_c == 0
  assert cigar_c == 1
  assert d == -11


def analyze_read_test8():
  """analyze read: check read with insertion in the middle"""
  r = pysam.AlignedSegment()
  r.qname = '1|1|0|0|21|100|30=20I50=|1|201|100|100='
  r.is_paired = True
  r.is_read1 = True

  r.cigarstring = '30M20I50M'
  r.reference_id = 0
  r.pos = 21
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d, long_insert, long_insert_m = creed.analyze_read(r)

  assert chrom_c == 1
  assert pos_c == 1
  assert cigar_c == 1
  assert d == 0


def analyze_read_test9():
  """analyze read: check read with insertion in the middle and position across the breakpoint"""
  r = pysam.AlignedSegment()
  r.qname = '1|1|0|0|21|100|30=20I50=|1|201|100|100='
  r.is_paired = True
  r.is_read1 = True

  r.cigarstring = '30M20I50M'
  r.reference_id = 0
  r.pos = 51
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d, long_insert, long_insert_m = creed.analyze_read(r)

  assert chrom_c == 1
  assert pos_c == 1
  assert cigar_c == 1
  assert d == 0

  # Read placed a little off, verify d is correct
  r.cigarstring = '30M20I50M'
  r.reference_id = 0
  r.pos = 61
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d, long_insert, long_insert_m = creed.analyze_read(r)

  assert chrom_c == 1
  assert pos_c == 0
  assert d == 10


def analyze_read_test10():
  """analyze read: check read with deletion in the middle and position across the breakpoint"""
  r = pysam.AlignedSegment()
  r.qname = '1|1|0|0|21|100|30=20D70=|1|201|100|100='
  r.is_paired = True
  r.is_read1 = True

  r.cigarstring = '30M20D70M'
  r.reference_id = 0
  r.pos = 21
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d, long_insert, long_insert_m = creed.analyze_read(r)

  assert chrom_c == 1
  assert pos_c == 1
  assert cigar_c == 1
  assert d == 0

  r.cigarstring = '30M20D70M'
  r.reference_id = 0
  r.pos = 71
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d, long_insert, long_insert_m = creed.analyze_read(r)

  assert chrom_c == 1
  assert pos_c == 1
  assert cigar_c == 1
  assert d == 0

  r.cigarstring = '30M20D70M'
  r.reference_id = 0
  r.pos = 81
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d, long_insert, long_insert_m = creed.analyze_read(r)

  assert chrom_c == 1
  assert pos_c == 0
  assert cigar_c == 1
  assert d == 10


def analyze_read_test11():
  """analyze read: long deletion in middle pos across breakpoint"""
  r = pysam.AlignedSegment()
  r.qname = '716964|5|1|0|4465332|100|33=250D67=|1|4465985|100|100='
  r.is_paired = True
  r.is_read1 = True

  r.cigarstring = '33S53M14S'
  r.reference_id = 4
  r.pos = 4465615
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d, long_insert, long_insert_m = creed.analyze_read(r)

  assert chrom_c == 1
  assert pos_c == 1
  assert cigar_c == 0
  assert d == 0
