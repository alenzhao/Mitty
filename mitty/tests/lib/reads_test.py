import numpy as np

import mitty.lib.variants as vr
import mitty.lib.reads as reads
from nose.tools import assert_sequence_equal


def expand_seq_test1():
  """Sequence expand no variants"""
  ref_seq = 'ACTGACTGACTGACT'
  ml = vr.VariantList()
  chrom = []
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 0)

  assert ref_seq == alt_seq
  assert_sequence_equal(beacons[1:-1], [])


def expand_seq_test2():
  """Sequence expand one variant of each type"""
  #          012345678901234
  ref_seq = 'ACTGACTGACTGACT'

  pos = [3]
  stop = [4]
  ref = ['G']
  alt = ['T']
  p = [0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 1)]
  #        012345678901234
  #ref     ACTGACTGACTGACT
  m_alt = 'ACTTACTGACTGACT'
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 0)
  assert ref_seq == alt_seq, alt_seq
  assert_sequence_equal(beacons[1:-1], [], str(beacons))

  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 1)
  assert alt_seq == m_alt, alt_seq
  assert_sequence_equal(beacons[1:-1].tolist(), [(3, 3, 0)])

  pos = [3]
  stop = [4]
  ref = ['G']
  alt = ['GAAA']
  p = [0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 0)]

  #        0123   45678901234
  #ref     ACTG   ACTGACTGACT
  m_alt = 'ACTGAAAACTGACTGACT'
  #        012345678901234567
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 0)
  assert alt_seq == m_alt, alt_seq
  assert_sequence_equal(beacons[1:-1].tolist(), [(4, 4, 3)])

  pos = [3]
  stop = [7]
  ref = ['GACT']
  alt = ['G']
  p = [0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 2)]

  #        012345678901234
  #            xxx
  #ref     ACTGACTGACTGACT
  m_alt = 'ACTGGACTGACT'
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 0)
  assert alt_seq == m_alt, alt_seq
  assert_sequence_equal(beacons[1:-1].tolist(), [(4, 4, -3)])


def expand_seq_test3():
  """Sequence expand multiple variants, different combinations"""
  #          123456789012
  ref_seq = 'ACTGACTGACTG'

  pos = [3, 5]
  stop = [4, 8]
  ref = ['G', 'CTG']
  alt = ['GAT', 'C']
  p = [0.1, 0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 1), (1, 1)]
  #        1234  5678   9012
  #ref     ACTG  ACTG   ACTG
  m_alt = 'ACTGATAC' + 'ACTG'
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 1)
  assert alt_seq == m_alt, alt_seq
  assert_sequence_equal(beacons[1:-1].tolist(), [(4, 4, 2), (6, 8, -2)], beacons)

  chrom = [(0, 0), (1, 1)]
  #        12345678   9012
  #ref     ACTGACTG   ACTG
  m_alt = 'ACTGAC' + 'ACTG'
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 1)
  assert alt_seq == m_alt, alt_seq
  assert_sequence_equal(beacons[1:-1].tolist(), [(6, 6, -2)], beacons)

  pos = [1, 5, 7]
  stop = [4, 6, 8]
  ref = ['CTG', 'C', 'G']
  alt = ['C', 'T', 'GAT']
  p = [0.1, 0.1, 0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 2), (1, 0), (2, 2)]
  #        12345678  9012
  #ref     ACTGACTG  ACTG
  m_alt = 'AC''ATTGATACTG'
  #        12  3456789012
  #         ^   ^ ^
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 0)
  assert alt_seq == m_alt, alt_seq
  assert_sequence_equal(beacons[1:-1].tolist(), [(2, 2, -2), (5, 3, 0), (8, 6, 2)])


def cigar_test1():
  """Rolling cigars: No variants"""
  #          012345678901
  ref_seq = 'ACTGACTGACTG'
  ml = vr.VariantList()
  chrom = []
  alt_seq, variant_waypoints = reads.expand_sequence(ref_seq, ml, chrom, 0)
  read_list = np.rec.fromrecords([(n, 3) for n in range(0, 8)], names=['start_a', 'read_len'])
  cigars = reads.roll_cigars(variant_waypoints, read_list)
  for cigar in cigars:
    assert cigar == '3M', cigar


def cigar_test2():
  """Rolling cigars: SNP"""
  #          012345678901
  ref_seq = 'ACTGACTGACTG'
  #             ^
  pos = [3]
  stop = [4]
  ref = ['G']
  alt = ['T']
  p = [0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 2)]

  #      012345678901234
  #ref   ACTGACTGACTGACT
  #alt   ACTTACTGACTGACT
  #         ^
  alt_seq, variant_waypoints = reads.expand_sequence(ref_seq, ml, chrom, 0)
  read_list = np.rec.fromrecords([(n, 3) for n in range(6)], names=['start_a', 'read_len'])
  cigars = reads.roll_cigars(variant_waypoints, read_list)
  assert cigars[0] == '3M', cigars[0]
  assert cigars[1] == '2M1X', cigars[1]
  assert cigars[2] == '1M1X1M', cigars[2]
  assert cigars[3] == '1X2M', cigars[3]
  assert cigars[4] == '3M', cigars[4]
  assert cigars[5] == '3M', cigars[5]


def cigar_test3():
  """Rolling cigars: INS"""
  ref_seq = 'ACTGACTGACTG'
  pos = [3]
  stop = [4]
  ref = ['G']
  alt = ['GA']
  p = [0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 2)]

  #        0123 45678901234
  #ref     ACTG ACTGACTGACT
  #alt     ACTGAACTGACTGACT
  #        0123456789012345
  alt_seq, variant_waypoints = reads.expand_sequence(ref_seq, ml, chrom, 0)
  read_list = np.rec.fromrecords([(n, 3) for n in range(6)], names=['start_a', 'read_len'])
  cigars = reads.roll_cigars(variant_waypoints, read_list)
  assert cigars[0] == '3M', cigars[0]
  assert cigars[1] == '3M', cigars[1]
  assert cigars[2] == '2M1S', cigars[2]
  assert cigars[3] == '1M1I1M', cigars[3]
  assert cigars[4] == '1S2M', cigars[4]
  assert cigars[5] == '3M', cigars[5]


def cigar_test4():
  """Rolling cigars: read completely inside INS"""
  ref_seq = 'ACTGACTGACTG'
  pos = [3]
  stop = [4]
  ref = ['G']
  alt = ['GAAAA']
  p = [0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 2)]

  #        0123    45678901234
  #ref     ACTG    ACTGACTGACT
  #alt     ACTGAAAAACTGACTGACT
  #        0123456789012345678
  alt_seq, variant_waypoints = reads.expand_sequence(ref_seq, ml, chrom, 0)
  read_list = np.rec.fromrecords([(n, 3) for n in range(7)], names=['start_a', 'read_len'])
  cigars = reads.roll_cigars(variant_waypoints, read_list)
  assert cigars[0] == '3M', cigars[0]
  assert cigars[1] == '3M', cigars[1]
  assert cigars[2] == '2M1S', cigars[2]
  assert cigars[3] == '1M2S', cigars[3]
  assert cigars[4] == '3S', cigars[4]
  assert cigars[5] == '3S', cigars[5]
  assert cigars[6] == '2S1M', cigars[6]