"""Tests for the generation of reads by applying variants to a reference"""
import os

from nose.tools import assert_raises, assert_sequence_equal
import numpy as np
from numpy.testing import assert_array_equal

import mitty.test
import mitty.simulation.readgen as rgen
import mitty.lib.vcfio as vio


def test_v_type():
  """Read gen: Flag complex variants in VCF"""
  df = vio.read_sample_from_vcf(os.path.join(mitty.test.example_data_dir, 'flawed-tiny.vcf'), 'g0_s0')
  for _, v in df.iterrows():
    assert_raises(RuntimeError, rgen.get_v_type, v)


def load_data():
  """A set of reference and variants for testing read generation."""
  seq = open(os.path.join(mitty.test.example_data_dir, 'tiny.fasta')).readlines()[1]
  df = vio.read_sample_from_vcf(os.path.join(mitty.test.example_data_dir, 'tiny.vcf'), 'g0_s0')
  return seq, df


def test_snp_expansion1():
  """Read gen: SNP expansion basic"""
  ref_seq, df = load_data()
  v = df.ix[0]
  vt, vl = rgen.get_v_type(v)
  assert vt == 'X', vt
  assert vl == 0, vl


def test_snp_expansion2():
  """Read gen: SNP expansion: No M section"""
  ref_seq, df = load_data()
  v = df.ix[0]
  samp_pos = 1
  ref_pos = 5
  vt, vl = rgen.get_v_type(v)
  nodes, n_samp_pos, n_ref_pos = rgen.snp(ref_seq, samp_pos, ref_pos, v, vl)

  assert len(nodes) == 1, nodes
  assert nodes[0] == (1, 5, 'X', 1, 'T', 0), nodes[0]
  assert n_samp_pos == samp_pos + 1
  assert n_ref_pos == ref_pos + 1


def test_snp_expansion3():
  """Read gen: SNP expansion: M section"""
  ref_seq, df = load_data()
  v = df.ix[0]
  samp_pos = 1
  ref_pos = 1
  vt, vl = rgen.get_v_type(v)
  nodes, samp_pos, ref_pos = rgen.snp(ref_seq, samp_pos, ref_pos, v, vl)

  assert len(nodes) == 2, nodes
  assert nodes[0] == (1, 1, '=', 4, 'ATGA', None), nodes
  assert nodes[1] == (5, 5, 'X', 1, 'T', 0), nodes


def test_ins_expansion1():
  """Read gen: INS expansion basic"""
  ref_seq, df = load_data()
  v = df.ix[1]
  vt, vl = rgen.get_v_type(v)
  assert vt == 'I', vt
  assert vl == 3, vl


def test_ins_expansion2():
  """Read gen: INS expansion: No M section"""
  ref_seq, df = load_data()
  v = df.ix[1]
  samp_pos = 9
  ref_pos = 9
  vt, vl = rgen.get_v_type(v)
  nodes, n_samp_pos, n_ref_pos = rgen.insertion(ref_seq, samp_pos, ref_pos, v, vl)

  assert len(nodes) == 1, nodes
  assert nodes[0] == (9, 9, 'I', 3, 'TTT', 3), nodes[0]
  assert n_samp_pos == samp_pos + 3
  assert n_ref_pos == ref_pos


def test_ins_expansion3():
  """Read gen: INS expansion: M section"""
  ref_seq, df = load_data()
  v = df.ix[1]
  samp_pos = 6
  ref_pos = 6
  vt, vl = rgen.get_v_type(v)
  nodes, samp_pos, ref_pos = rgen.insertion(ref_seq, samp_pos, ref_pos, v, vl)

  assert len(nodes) == 2, nodes
  assert nodes[0] == (6, 6, '=', 3, 'GTA', None), nodes[0]
  assert nodes[1] == (9, 9, 'I', 3, 'TTT', 3), nodes[1]


def test_del_expansion1():
  """Read gen: DEL expansion basic"""
  ref_seq, df = load_data()
  v = df.ix[2]
  vt, vl = rgen.get_v_type(v)
  assert vt == 'D', vt
  assert vl == 2, vl


def test_del_expansion2():
  """Read gen: DEL expansion: No M section"""
  ref_seq, df = load_data()
  v = df.ix[2]
  samp_pos = 12
  ref_pos = 12
  vt, vl = rgen.get_v_type(v)
  nodes, n_samp_pos, n_ref_pos = rgen.deletion(ref_seq, samp_pos, ref_pos, v, vl)

  assert len(nodes) == 1, nodes
  assert nodes[0] == (11, 14, 'D', 2, '', -2), nodes[0]
  assert n_samp_pos == samp_pos, n_samp_pos
  assert n_ref_pos == ref_pos + 2, n_ref_pos


def test_del_expansion3():
  """Read gen: DEL expansion: M section"""
  ref_seq, df = load_data()
  v = df.ix[2]
  samp_pos = 12
  ref_pos = 9
  vt, vl = rgen.get_v_type(v)
  nodes, samp_pos, ref_pos = rgen.deletion(ref_seq, samp_pos, ref_pos, v, vl)

  assert len(nodes) == 2, nodes
  assert nodes[0] == (12, 9, '=', 3, 'TCC', None), nodes
  assert nodes[1] == (14, 14, 'D', 2, '', -2), nodes


def test_expand_sequence():
  """Read gen: Sequence expansion"""
  ref_seq, df = load_data()
  nodes = rgen.create_node_list(ref_seq, ref_start_pos=1, chrom_copy=0b01, vcf=df)

  assert len(nodes) == 9
  #                   ps pr op l seq
  assert nodes[0] == (1, 1, '=', 4, 'ATGA', None)
  assert nodes[1] == (5, 5, 'X', 1, 'T', 0)
  assert nodes[2] == (6, 6, '=', 3, 'GTA', None)
  assert nodes[3] == (9, 9, 'I', 3, 'TTT', 3)
  assert nodes[4] == (12, 9, '=', 3, 'TCC', None)
  assert nodes[5] == (14, 14, 'D', 2, '', -2)
  assert nodes[6] == (15, 14, '=', 7, 'GGAGGCG', None)
  assert nodes[7] == (21, 23, 'D', 2, '', -2)
  assert nodes[8] == (22, 23, '=', 3, 'ACC', None)


def test_get_begin_end_nodes():
  """Read gen: find start and stop nodes for reads"""
  ref_seq, df = load_data()
  nodes = rgen.create_node_list(ref_seq, ref_start_pos=1, chrom_copy=0b01, vcf=df)

  pl = np.arange(1, 16, dtype=int)
  ll = 10
  n0, n1 = rgen.get_begin_end_nodes(pl, ll, nodes)
  assert_array_equal(n0, np.array([0, 0, 0, 0, 1, 2 , 2, 2, 3, 3, 3, 4, 4, 4, 6]))
  assert_array_equal(n1, np.array([3, 3, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 8, 8, 8]))


def test_read_gen1():
  """Read gen: Read pos, cigar, v_list and seq (cpy 1)"""
  ref_seq, df = load_data()
  nodes = rgen.create_node_list(ref_seq, ref_start_pos=1, chrom_copy=0b01, vcf=df)

  assert rgen.generate_read(1, 10, 0, 3, nodes) == (1, '4=1X3=2I', [0, 3], 'ATGATGTATT')
  assert rgen.generate_read(2, 10, 0, 3, nodes) == (2, '3=1X3=3I', [0, 3], 'TGATGTATTT')
  assert rgen.generate_read(3, 10, 0, 4, nodes) == (3, '2=1X3=3I1=', [0, 3], 'GATGTATTTT')
  assert rgen.generate_read(4, 10, 0, 4, nodes) == (4, '1=1X3=3I2=', [0, 3], 'ATGTATTTTC')
  assert rgen.generate_read(5, 10, 1, 4, nodes) == (5, '1X3=3I3=', [0, 3], 'TGTATTTTCC')
  assert rgen.generate_read(6, 10, 2, 6, nodes) == (6, '3=3I3=2D1=', [3, -2], 'GTATTTTCCG')
  assert rgen.generate_read(7, 10, 2, 6, nodes) == (7, '2=3I3=2D2=', [3, -2], 'TATTTTCCGG')
  assert rgen.generate_read(8, 10, 2, 6, nodes) == (8, '1=3I3=2D3=', [3, -2], 'ATTTTCCGGA')
  assert rgen.generate_read(9, 10, 3, 6, nodes) == (9, '3I3=2D4=', [3, -2], 'TTTTCCGGAG')
  assert rgen.generate_read(10, 10, 3, 6, nodes) == (9, '2I3=2D5=', [3, -2], 'TTTCCGGAGG')
  assert rgen.generate_read(11, 10, 3, 6, nodes) == (9, '1I3=2D6=', [3, -2], 'TTCCGGAGGC')
  assert rgen.generate_read(12, 10, 4, 6, nodes) == (9, '3=2D7=', [-2], 'TCCGGAGGCG')
  assert rgen.generate_read(13, 10, 4, 8, nodes) == (10, '2=2D7=2D1=', [-2, -2], 'CCGGAGGCGA')
  assert rgen.generate_read(14, 10, 4, 8, nodes) == (11, '1=2D7=2D2=', [-2, -2], 'CGGAGGCGAC')
  assert rgen.generate_read(15, 10, 6, 8, nodes) == (14, '7=2D3=', [-2], 'GGAGGCGACC')


def test_read_gen2():
  """Read gen: Read pos, cigar, v_list and seq (cpy 2)"""
  ref_seq, df = load_data()
  nodes = rgen.create_node_list(ref_seq, ref_start_pos=1, chrom_copy=0b10, vcf=df)

  assert rgen.generate_read(1, 10, 0, 0, nodes) == (1, '10=', [], 'ATGACGTATC')
  assert rgen.generate_read(2, 10, 0, 0, nodes) == (2, '10=', [], 'TGACGTATCC')
  assert rgen.generate_read(4, 10, 0, 0, nodes) == (4, '10=', [], 'ACGTATCCAA')
  assert rgen.generate_read(5, 10, 0, 1, nodes) == (5, '9=1X', [0], 'CGTATCCAAT')
  assert rgen.generate_read(6, 10, 0, 2, nodes) == (6, '8=1X1=', [0], 'GTATCCAATG')
  # TODO. Add more test cases? The core test - whether we take the right variant - is already tested by these cases


def test_read_gen_ins():
  """Read gen: Reads from inside insertion"""
  ref_seq, df = load_data()
  nodes = rgen.create_node_list(ref_seq, ref_start_pos=1, chrom_copy=0b01, vcf=df)

  assert rgen.generate_read(9, 2, 3, 3, nodes) == (8, '>0:2I', [3], 'TT')