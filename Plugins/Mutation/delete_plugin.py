"""This is the stock deletion generator. Please see Readme for details on how to write variant generator plugins for
mutate.py.

The plugin uses two random number generators. The first creates poisson distributed numbers to locate the deletions.
The other creates poisson distributed numbers to determine deletion lengths. This is equivalent to modeling deletion
termination as a bernoulli process.

Note: This never generates a deletion at the first base of a sequence.

Example parameter snippet

    "mydelete": {
        "model": "delete",
        "start_dels_frac": 0.7,
        "stop_dels_frac":  0.9,
        "phet": 0,
        "p_del": 0.01,
        "lam_del": 10,
        "het_rng_seed": 3,
        "strand_rng_seed": 4,
        "del_loc_rng_seed": 0,
        "del_len_rng_seed": 1
    }
"""
import numpy
import logging
het_type = ['0/1', '1/0']  # The two types of het

logger = logging.getLogger(__name__)


def variant(ref_seq=None, ref_seq_len=0,
            phet=0.0, p_del=0.01, lam_del=5,
            start_dels_frac=0.0, stop_dels_frac=1.0,
            het_rng_seed=1,
            strand_rng_seed=4,
            del_loc_rng_seed=1,
            del_len_rng_seed=1,
            block_size=10000, **kwargs):
  """A generator which returns a variant when asked for. This is the stock SNP generator and returns snp locations in
  a poisson distributed fashion.
  Inputs:
    ref_seq              - The reference sequence
    ref_seq_len          - length of whole sequence (needed to compute start and stop)
    phet                 - probability of having heterozygous mutation
    p_del                - probability of deletes
    lam_del              - mean length of poisson distributed delete lengths
    start_dels_frac      - start generating dels from here (0.0, 1.0)
    stop_dels_frac       - stop generating dels after this (0.0, 1.0) stop_snps_frac > start_snps_frac
    het_rng_seed         - rng used to decide if genotype is heterozygous or not (0/1 or 1/0  vs 1/1)
    del_loc_rng_seed     - SNP locator rng numpy.random.RandomState(seed)
    del_len_rng_seed     - rng used to determine length of delete
    kwargs               - absorbs any other parameters it does not use

  Outputs:
    variant              - (POS, REF, ALT, GT, skip, list(footprints))


  Test with 'N's. No deletions should straddle a region with N
  >>> args = {'p_del': .1, 'lam_del': 3, 'del_loc_rng_seed': 1, 'del_len_rng_seed': 2}; \
  ref_seq='ACGTACGTANGTACGTACGTACGTACGTACGTACGTACGTACNTACGTACGTACGTACGT'; ref_seq_len = len(ref_seq); \
  gen = variant(ref_seq=ref_seq, ref_seq_len=ref_seq_len, block_size=100, **args);
  >>> for n in range(10): print next(gen,None)
  (15, 'TACG', 'T', '1/1', 20, None)
  (22, 'GTA', 'G', '1/1', 26, None)
  (31, 'TACG', 'T', '1/1', 36, None)
  (49, 'CGT', 'C', '1/1', 53, None)
  (55, 'TA', 'T', '1/1', 58, None)
  None
  None
  None
  None
  None

  Test with one block
  >>> args = {'p_del': .1, 'lam_del': 3, 'del_loc_rng_seed': 1, 'del_len_rng_seed': 2}; \
  ref_seq='ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'; ref_seq_len = len(ref_seq); \
  gen = variant(ref_seq=ref_seq, ref_seq_len=ref_seq_len, block_size=100, **args);
  >>> for n in range(10): print next(gen,None)
  (9, 'CG', 'C', '1/1', 12, None)
  (15, 'TACG', 'T', '1/1', 20, None)
  (22, 'GTA', 'G', '1/1', 26, None)
  (31, 'TACG', 'T', '1/1', 36, None)
  (40, 'ACGTA', 'A', '1/1', 46, None)
  (49, 'CGT', 'C', '1/1', 53, None)
  (55, 'TA', 'T', '1/1', 58, None)
  None
  None
  None

  Test with multiple blocks - should be same answer even though we have changed the block size
  >>> gen = variant(ref_seq=ref_seq, ref_seq_len=ref_seq_len, block_size=1, **args);
  >>> for n in range(10): print next(gen,None)
  (9, 'CG', 'C', '1/1', 12, None)
  (15, 'TACG', 'T', '1/1', 20, None)
  (22, 'GTA', 'G', '1/1', 26, None)
  (31, 'TACG', 'T', '1/1', 36, None)
  (40, 'ACGTA', 'A', '1/1', 46, None)
  (49, 'CGT', 'C', '1/1', 53, None)
  (55, 'TA', 'T', '1/1', 58, None)
  None
  None
  None

  Test heterozygous deletes
  >>> args = {'phet': 0.5, 'p_del': .1, 'lam_del': 3, 'het_rng_seed': 3, 'strand_rng_seed': 5, 'del_loc_rng_seed': 1, 'del_len_rng_seed': 2}; \
  ref_seq='ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'; ref_seq_len = len(ref_seq); \
  gen = variant(ref_seq=ref_seq, ref_seq_len=ref_seq_len, block_size=100, **args);
  >>> for n in range(10): print next(gen,None)
  (9, 'CG', 'C', '1/1', 12, None)
  (15, 'TACG', 'T', '1/1', 20, None)
  (22, 'GTA', 'G', '1/0', 26, None)
  (31, 'TACG', 'T', '1/1', 36, None)
  (40, 'ACGTA', 'A', '1/1', 46, None)
  (49, 'CGT', 'C', '1/1', 53, None)
  (55, 'TA', 'T', '0/1', 58, None)
  None
  None
  None
  """
  def get_locs_and_lens():  # Simply a convenience.
    het_or_not = het_rng.rand(block_size)
    strand_no = strand_rng.randint(2, size=block_size)
    loc_diff = numpy.maximum(del_loc_rng.poisson(lam=1.0 / p_del, size=block_size), 1)
    del_lens = del_len_rng.poisson(lam=lam_del, size=block_size)
    locs = numpy.cumsum(loc_diff) + start_offset
    return het_or_not, strand_no, locs, del_lens, locs[-1]

  het_rng = numpy.random.RandomState(seed=het_rng_seed)
  strand_rng = numpy.random.RandomState(seed=strand_rng_seed)
  del_loc_rng = numpy.random.RandomState(seed=del_loc_rng_seed)
  del_len_rng = numpy.random.RandomState(seed=del_len_rng_seed)

  start_offset = int(ref_seq_len * start_dels_frac)
  del_end = int(ref_seq_len * stop_dels_frac)

  het_or_not, strand_no, locs, del_lens, start_offset = get_locs_and_lens()
  internal_cntr = 0
  while locs[internal_cntr] < del_end:
    vl = locs[internal_cntr]
    dl = del_lens[internal_cntr]
    ref = ref_seq[vl:vl+dl+1]
    if 'N' in ref:
      alt = None  # Very conservative - we only do deletions in completely known regions
    else:
      alt = ref[0]
      gt = '1/1' if het_or_not[internal_cntr] > phet else het_type[strand_no[internal_cntr]]

    internal_cntr += 1
    if internal_cntr == locs.size:
      het_or_not, strand_no, locs, del_lens, start_offset = get_locs_and_lens()
      internal_cntr = 0

    if alt is not None:
      yield (vl, ref, alt, gt, vl + dl + 2, None)  # POS, REF, ALT, skipto, list(footprints)
                                               # footprints, in this case, is None, since we simply skip forward
      # We have vl dl + 2 because we want 1 base buffer between variants, even SNPs (see Readme)

if __name__ == "__main__":
  import doctest
  doctest.testmod()
