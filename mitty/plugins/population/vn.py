"""
DESCRIPTION:

The vn model was designed to generate one sample genome (S) and a series of graphs (G0, G1, ...)

S has a fraction of variants not found in G0, G1 ... (novel variants)
and the rest are found in G0 (and G1, G2 ...) (known variants).

In set notation:

  S - G0 = S - G1 = S - G2 ... = novel variants
  S & G0 = S & G1 ... = known variants

  G0, G1, ... are designed like babushka dolls such that G0 is a subset of G1
  G1 a subset of G2 and so on

In set notation:

  G0 < G1 < G2 ... (< standing for subset)

WARNING:
Though, technically, it is possible to use G0, G1 ... as 'samples' to generate reads
this will be incorrect to do: G0, G1 may contain 'illegal' combinations of variations
such as a deletion with a SNP in the middle of it because they represent population graphs.
Such 'illegal' patterns lead to undefined behavior in the read generator.
The best you can hope for is an explicit error message in case of such use.

"""
import numpy as np
from collections import OrderedDict

__example_param_text = """
{
  "vn": {
    "p_S": 0.104,
    "p_G": [0.1, 0.2, 0.3],
  }
}
"""

__param_hints__ = """
PARAMETER SETTING HINTS:

The plugin ignores any site frequency spectrum model you set.
The density of the master list of variants derives directly from the variant model parameters.
p_S governs what fraction of the master list goes into the sample S
p_G values govern what fraction of the master list goes into the graphs G0, G1, ...
One constraint is the p_S > p_G[0] and p_S - p_G[0] governs the number of novel variants in
S that do not appear in G0, G1, G2, ...
"""

_description = __doc__ + '\nEXAMPLE PARAMETERS:\n' + __example_param_text + __param_hints__

#_example_params = json.loads(__example_param_text)
_example_params = eval(__example_param_text)


class Model:
  def __init__(self, p_S=0.104, p_G=[0.1, 0.2, 0.3]):
    """A population model that creates samples with more and more variants. Suitable for the aligner paper experiments

    :param p_S: probability value for S. Fraction of novel variants = p_S - p_G[0]
    :param p_G: probability values for g0, g1, g2 .... set
    """
    self.p_S, self.p_G = p_S, p_G

    assert 0 <= self.p_S <= 1.0, 'p_S needs to be >= 0 and <= 1.0'
    assert self.p_S > self.p_G[0], 'p_S needs to be > p_G[0]'
    for n in range(len(self.p_G) - 1):
      assert self.p_G[n] < self.p_G[n + 1], 'p_G needs to be in ascending order'
      assert 0 <= self.p_G[n] <= 1.0, 'p_G needs to be >= 0 and <= 1.0'

  def samples(self, chrom_no=None, ml=None, rng_seed=1, **kwargs):
    """This returns an iterator

    :param chrom_no:  number of the chromosome being considered [1,2,3 ...]  (ignored here)
    :param ml:        VariantList. master list of variants as created by genomes program
    :param rng_seed:  seed for random number generators
    :return: A generator returning (generation no, serial_no, chromosome, % samples done) for each sample in population

    Algorithm:

    Generate random numbers r same size as master list
    r < p_S => S
    r < p_G[0] => G0
    p_G[0] < r < p_S = S - G0 (The novels). Set these to > 1 and take them out of circulation
    r < p_G[1] => G1 and so on.

    For p_S, pick genotype using uniform random numbers such that p([0|1]) = p([1|0]) = p([1|1]) = 1/3
    """

    rng = np.random.RandomState(rng_seed)
    r = rng.rand(ml.variants.shape[0])
    for n in range(len(self.p_G) + 1):
      if n == 0:
        sample_name = 'S'
        v_idx = (r < self.p_S).nonzero()[0]
        gt = rng.randint(0, 3, v_idx.size)
        chrom = ml.zip_up_chromosome(v_idx[(gt == 0) | (gt == 2)], v_idx[(gt == 1) | (gt == 2)])
      elif n == 1:
        sample_name = 'G0'
        v_idx = (r < self.p_G[0]).nonzero()[0]
        chrom = np.empty(shape=(v_idx.size,), dtype=[('index', 'i4'), ('gt', 'i1')])
        chrom['index'] = v_idx
        chrom['gt'] = 2  # Graph genotype is
        r[(self.p_G[0] <= r) & (r < self.p_S)] = 1.1  # Take S - G0 (The novels) out of circulation
      else:
        sample_name = 'G{}'.format(n - 1)
        v_idx = (r < self.p_G[n - 1]).nonzero()[0]
        chrom = np.empty(shape=(v_idx.size,), dtype=[('index', 'i4'), ('gt', 'i1')])
        chrom['index'] = v_idx
        chrom['gt'] = 2

      yield sample_name, chrom, float(n + 1) / self.get_sample_count_estimate()

  def get_sample_count_estimate(self):
    """Give us an as exact as possible estimate of how many samples we will produce"""
    return 1 + len(self.p_G)

  def inspect(self, pop):
    """Given a population created using this plugin in return a text string with a description of what this population
    is about."""
    chr_list = pop.get_chromosome_list()
    v_list = ['S'] + ['G{:d}'.format(n) for n in range(len(self.p_G))]

    counts = OrderedDict(
      [(k, []) for k in v_list] +
      [(k, []) for i in range(1, len(v_list)) for k in [s.format(v_list[0], v_list[i]) for s in ['{} - {}', '{} & {}']]]
    )

    for n_chr, chr in enumerate(chr_list):
      v_idx = [set(pop.get_sample_variant_index_for_chromosome(chr, v)['index']) for v in v_list]
      for i, v in enumerate(v_list):
        counts[v] += [len(v_idx[i])]
      for i in range(1, len(v_list)):
        counts['{} - {}'.format(v_list[0], v_list[i])] += [len(v_idx[0] - v_idx[i])]
        counts['{} & {}'.format(v_list[0], v_list[i])] += [len(v_idx[0] & v_idx[i])]

    rep_str = __doc__
    rep_str += "VARIANT COUNTS FOR THIS DB:\n\n"
    row_format = "    " + "{:>12} " * (len(counts) + 1)
    sep = "    " + "-" * 13 * (len(counts) + 1) + '\n'
    rep_str += row_format.format("chrom", *["|{}|".format(k) for k in counts.keys()]) + '\n'
    rep_str += sep
    for n, chrom in enumerate(chr_list):
      rep_str += row_format.format(chrom, *[cnts[n] for cnts in counts.values()]) + '\n'
    rep_str += sep
    rep_str += row_format.format("Total", *[sum(cnts) for cnts in counts.values()]) + '\n'

    return rep_str
