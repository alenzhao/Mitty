"""This is the stock deletion generator. The length of deletions follows a geometric distribution as would be expected
from a poisson point process governing deletion termination"""
import itertools

import numpy as np

import mitty.lib
import mitty.lib.util as mutil
from mitty.plugins.variants import scale_probability_and_validate

import logging
logger = logging.getLogger(__name__)


__example_param_text__ = """
{
  "p": 0.01,           # probability that the deletion will happen at any given base
  "p_end": 0.1,        # probability governing length of deletion
  "del_len_max": 1000  # upper bound on deletion lengths
}
"""

_description = """
This is the stock delete plugin. A typical parameter set resembles
""" + __example_param_text__

_example_params = eval(__example_param_text__)


class Model:
  def __init__(self, p=0.01, p_end=0.1, del_len_max=1000, **kwargs):
    assert 0 <= p <= 1.0, "Probability out of range"
    assert 0 <= p_end <= 1.0, "Probability out of range"
    assert 0 < del_len_max, "Maximum deletion length needs to be greater than 0"
    self.p, self.p_end, self.del_len_max = p, p_end, del_len_max

  def get_variants(self, ref, chrom, p, f, seed=1):
    """This function is called by the simulator to obtain variants.

    :param ref: reference sequence as a string
    :param chrom: chromosome number (1,2,3,4...)
    :param p: array/list of probability values
    :param f: array/list of frequency values
    :param seed: seed for the random number generators
    :return: 5 arrays/lists/iterables all of the same length
              pos   - position of SNPs
              stop  - stop locations, (pos + 1 for SNPs)
              ref   - reference base,
              alt   - alt base,
              p     - probability value for this variant. These are uniformly distributed random values
    """
    assert 0 < seed < mitty.lib.SEED_MAX
    logger.debug('Master seed: {:d}'.format(seed))

    base_loc_rng, del_len_rng = mutil.initialize_rngs(seed, 2)

    p_eff = scale_probability_and_validate(self.p, p, f)
    del_locs = mutil.place_poisson_seq(base_loc_rng, p_eff, 0, len(ref), ref)  #np.array([x for x in mutil.place_poisson(base_loc_rng, p_eff, 0, len(ref)) if ref[x] != 'N'], dtype='i4')
    del_lens = del_len_rng.geometric(p=p_eff, size=del_locs.shape[0])  #mutil.base_subs(ref, snp_locs, self.t_mat, base_t_rng)
    np.clip(del_lens, a_min=2, a_max=self.del_len_max, out=del_lens)  # Make sure our deletions are clipped at the level we want
    idx = ((del_locs + del_lens) < len(ref)).nonzero()[0]   # Get rid of any deletions that go past the sequence end
    del_locs = del_locs[idx]
    del_lens = del_lens[idx]
    # http://stackoverflow.com/questions/8081545/convert-list-of-tuples-to-multiple-lists-in-python
    idx, refs, alts = map(list, itertools.izip(*((n, ref[del_loc:del_loc + del_len], ref[del_loc]) for n, (del_loc, del_len) in enumerate(np.nditer([del_locs, del_lens])) if ref[del_loc + del_len - 1] != 'N')))
    # This gets rid of any deletions that stretch into the 'N' regions of a sequence
    return del_locs[idx], del_locs[idx] + del_lens[idx], refs, alts, del_lens[idx] / float(del_lens[idx].max())


def test():
  """Basic test"""
  ref_seq = 'ACTGACTGACTGACTGACTGACTGACTGACTGACTG'
  m = Model(p=0.1)
  pos, stop, ref, alt, p = m.get_variants(ref_seq, 1, np.array([0.2]), np.array([1.0]), seed=10)
  for p, r in zip(pos, alt):
    assert r == ref_seq[p]


if __name__ == "__main__":
  print _description