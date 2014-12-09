"""This is the stock SNP plugin. It uses four independent RNGs to locate SNPs along a reference sequence, assign each
SNP a zygosity and assign an ALT base.
"""
from mitty.lib.variation import new_variation
import mitty.lib.util as mutil
import logging
logger = logging.getLogger(__name__)

__example_param_text = """
{
  "chromosome": [4],      # List of chromosomes to apply the variant to
  "p": 0.01,              # probability that the SNP will happen at any given base
  "phet": 0.5,            # probability that the variant will be heterozygous
}
"""

_description = """
This is the stock SNP plugin. A typical parameter set resembles
""" + __example_param_text

_example_params = eval(__example_param_text)


# This is a substitution base substitution table
base_sub_dict = {
  'A': ['T', 'C', 'G'],  # 'A' -> 'TCG'
  'C': ['G', 'A', 'T'],  # 'C' -> 'GAT'
  'G': ['A', 'T', 'C'],  # 'G' -> 'ATC'
  'T': ['C', 'G', 'A'],  # 'T' -> 'CGA'
  'N': ['N', 'N', 'N']

}


def variant_generator(ref={},
                      chromosome=None,
                      phet=0.5,
                      p=0.01,
                      master_seed=1,
                      **kwargs):
  assert 0 <= p <= 1.0, "Probability out of range"
  assert 0 <= phet <= 1.0, "Probability out of range"
  logger.debug('Master seed: {:d}'.format(master_seed))
  base_loc_rng, base_sub_rng, het_rng, copy_rng = mutil.initialize_rngs(master_seed, 4)

  for chrom in chromosome:
    ref_chrom = ref[chrom]
    snp_locs = mutil.place_poisson(base_loc_rng, p, len(ref_chrom))
    het_type = mutil.zygosity(snp_locs.size, phet, het_rng, copy_rng)
    base_subs = base_sub_rng.randint(3, size=snp_locs.size)

    yield {chrom: [new_variation(pos + 1, pos + 2, ref_chrom[pos], base_sub_dict[ref_chrom[pos]][bs], het)
                   for pos, bs, het in zip(snp_locs, base_subs, het_type) if ref_chrom[pos] != 'N']}
    # +1 because VCF files are 1 indexed
    #alts will be 0 if ref is not one of ACTG


def test():
  from mitty.lib.variation import Variation
  """Basic test"""
  ref = {
    1: 'ACTGACTGACTG',
    2: 'ACTGACTGACTGACTGACTGACTG'
  }
  vg = variant_generator(ref, chromosome=[1, 2], p=1.0)
  for v_list in vg:
    assert isinstance(v_list.values()[0][0], Variation), v_list

if __name__ == "__main__":
  print _description