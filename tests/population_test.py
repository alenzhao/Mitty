from . import *
import numpy.testing
from mitty.variation import Variation
from mitty.population import *
import mitty.denovo as denovo
from nose.tools import raises


def chrom_crossover_test():
  c1 = [
      Variation(1, 2, 'C', 'CAA', HET2),
      Variation(3, 6, 'CAG', 'C', HET1),
      Variation(7, 8, 'G', 'T', HET2),
      Variation(9, 12, 'GTT', '', HET1),
      Variation(13, 16, 'GTT', 'TTG', HOMOZYGOUS),
      Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS)
  ]
  crossover_idx = [1, 1, 1, 0, 1, 0]
  correct_chrom = [
      Variation(1, 2, 'C', 'CAA', HET1),
      Variation(3, 6, 'CAG', 'C', HET2),
      Variation(7, 8, 'G', 'T', HET1),
      Variation(9, 12, 'GTT', '', HET1),
      Variation(13, 16, 'GTT', 'TTG', HOMOZYGOUS),
      Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS)
  ]
  c1_c = chrom_crossover(c1, crossover_idx)
  assert correct_chrom == c1_c, c1_c
  assert c1[0] == Variation(1, 2, 'C', 'CAA', HET2)  # We shouldn't be changing our original


def crossover_event_test():
  c1 = [
    Variation(1, 2, 'C', 'CAA', HET2),
    Variation(3, 6, 'CAG', 'C', HET1),
    Variation(7, 8, 'G', 'T', HET2)
  ]
  c2 = [
    Variation(2, 5, 'CAG', 'C', HET1),
    Variation(7, 8, 'G', 'T', HET2)
  ]
  g1 = {1: c1, 2: c2}
  crossover_idx = {1: [1, 0, 1], 2: [0, 1]}

  correct_g = {
    1: [
      Variation(1, 2, 'C', 'CAA', HET1),
      Variation(3, 6, 'CAG', 'C', HET1),
      Variation(7, 8, 'G', 'T', HET1)
    ],
    2: [
      Variation(2, 5, 'CAG', 'C', HET1),
      Variation(7, 8, 'G', 'T', HET1)
    ]
  }
  assert correct_g == crossover_event(g1, crossover_idx)


def pair_one_chrom_test():
  """Merging: test c1 longer than c2"""
  c1 = [
    Variation(1, 2, 'C', 'CAA', HET2),
    Variation(3, 6, 'CAG', 'C', HET1),
    Variation(17, 18, 'G', 'T', HET2),
    Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS),
  ]
  c2 = [
    Variation(7, 8, 'G', 'T', HET2),
  ]
  which_copy = (1, 0)
  correct_pairing = [
    Variation(1, 2, 'C', 'CAA', HET1),
    Variation(17, 18, 'G', 'T', HET1),
    Variation(23, 26, 'GTT', 'TTG', HET1),
  ]
  assert correct_pairing == pair_one_chrom(c1, c2, which_copy), pair_one_chrom(c1, c2, which_copy)


def pair_one_chrom_test2():
  """Merging: test c2 longer than c1"""
  c1 = [
    Variation(1, 2, 'C', 'CAA', HET2)
  ]
  c2 = [
    Variation(7, 8, 'G', 'T', HET1),
    Variation(17, 18, 'G', 'T', HET1),
  ]
  which_copy = (1, 0)
  correct_pairing = [
    Variation(1, 2, 'C', 'CAA', HET1),
    Variation(7, 8, 'G', 'T', HET2),
    Variation(17, 18, 'G', 'T', HET2),
  ]
  assert correct_pairing == pair_one_chrom(c1, c2, which_copy), pair_one_chrom(c1, c2, which_copy)


def pair_one_chrom_test3():
  """Merging: Comprehensive test"""
  c1 = [
    Variation(1, 2, 'C', 'CAA', HET2),  # Tests homozygosity
    Variation(3, 6, 'CAG', 'C', HET1),  # Tests both variants are not on the copies chosen
    Variation(17, 18, 'G', 'T', HET2),  # Test zipper (several variants should come from c2 before we get to this)
    Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS),  # Tests handling of homozygous variants
    Variation(29, 30, 'T', 'G', HOMOZYGOUS)  # Tests unequal var list lengths
  ]
  c2 = [
    Variation(1, 2, 'C', 'CAA', HET1),
    Variation(3, 6, 'CAG', 'C', HET2),
    Variation(7, 8, 'G', 'T', HET2),
    Variation(9, 12, 'GTT', '', HET1),
    Variation(15, 18, 'GTT', 'G', HET1),  # Test partial overlap on different copies (should not affect each other)
    Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS)
  ]
  which_copy = (1, 0)
  correct_pairing = [
    Variation(1, 2, 'C', 'CAA', HOMOZYGOUS),  # Tests homozygosity
    Variation(9, 12, 'GTT', '', HET2),
    Variation(15, 18, 'GTT', 'G', HET2),  # Test partial overlap on different copies (should not affect each other)
    Variation(17, 18, 'G', 'T', HET1),  # Test zipper (several variants should come from c2 before we get to this)
    Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS),  # Tests handling of homozygous variants
    Variation(29, 30, 'T', 'G', HET1)  # Tests unequal var list lengths
  ]
  assert correct_pairing == pair_one_chrom(c1, c2, which_copy)


def fertilize_one_test():
  c11 = [
    Variation(1, 2, 'C', 'CAA', HET1)
  ]
  c12 = [
    Variation(7, 8, 'G', 'T', HET2),
    Variation(17, 18, 'G', 'T', HET2),
  ]
  g1 = {1: c11, 2: c12}

  c21 = [
    Variation(1, 2, 'C', 'CAA', HET2)
  ]
  c22 = [
    Variation(7, 8, 'G', 'T', HET1),
    Variation(17, 18, 'G', 'T', HET1),
  ]
  g2 = {1: c21, 2: c22}
  which_copy = {1: (0, 1), 2: (1, 0)}

  correct_g3 = {
    1: [
      Variation(1, 2, 'C', 'CAA', HOMOZYGOUS)
    ],
    2: [
      Variation(7, 8, 'G', 'T', HOMOZYGOUS),
      Variation(17, 18, 'G', 'T', HOMOZYGOUS),
    ]
  }
  assert correct_g3 == fertilize_one(g1, g2, which_copy)


def place_crossovers_on_chrom_test():
  """Cross over location generator"""
  c1 = [
    Variation(1, 2, 'C', 'CAA', HET2),
    Variation(3, 6, 'CAG', 'C', HET1),
    Variation(17, 18, 'G', 'T', HET2),
    Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS),
    Variation(29, 30, 'T', 'G', HOMOZYGOUS)
  ]

  hot_spots = numpy.array([])  # No hotspots, no crossover
  numpy.testing.assert_array_equal(numpy.array([0, 0, 0, 0, 0]),  # Hot spot is narrow and over first variant
                                   place_crossovers_on_chrom(c1, hot_spots, numpy.random.RandomState(seed=1)))

  hot_spots = numpy.array([[1, 1, .5]])
  numpy.testing.assert_array_equal(numpy.array([1, 0, 0, 0, 0]),  # Hot spot is narrow and over first variant
                                   place_crossovers_on_chrom(c1, hot_spots, numpy.random.RandomState(seed=1)))

  hot_spots = numpy.array([[17, 1, .5]])
  numpy.testing.assert_array_equal(numpy.array([0, 0, 1, 0, 0]),  # Hot spot is narrow and over third variant
                                   place_crossovers_on_chrom(c1, hot_spots, numpy.random.RandomState(seed=1)))

  hot_spots = numpy.array([[1, 1, .5], [17, 1, .5]])
  numpy.testing.assert_array_equal(numpy.array([1, 0, 1, 0, 0]),  # Two narrow hot spots over first and third variants
                                   place_crossovers_on_chrom(c1, hot_spots, numpy.random.RandomState(seed=1)))

  hot_spots = numpy.array([[23, 1, 100], [17, 1, 100]])
  numpy.testing.assert_array_equal(numpy.array([1, 1, 1, 1, 1]),  # Super broad hotspots, covers all
                                   place_crossovers_on_chrom(c1, hot_spots, numpy.random.RandomState(seed=1)))

  # A proper test is actually to do something like this (with c1 set as above)
  # import pylab
  # hot_spots = numpy.array([[3, 1, 2], [23, 1, 5]])
  # rng = numpy.random.RandomState(seed=1)
  # data = numpy.concatenate([numpy.nonzero(place_crossovers_on_chrom(c1, hot_spots, rng))[0] for n in range(100)])
  # pylab.hist(data)
  # pylab.plot()
  # The denser your variant structure, the clearer is the sum of gaussian model


def spawn_test():
  """Spawn from parents (No denovo)."""
  g1 = {1: [Variation(1, 2, 'C', 'CAA', HET1)],
        2: [Variation(7, 8, 'G', 'T', HET2),
            Variation(17, 18, 'G', 'T', HET2)]}
  g2 = {1: [Variation(1, 2, 'C', 'CAA', HET2)],
        2: [Variation(7, 8, 'G', 'T', HET1),
            Variation(17, 18, 'G', 'T', HET1)]}

  # These hotspots ensure crossing over at these locii
  hot_spots = {1: numpy.array([[1, 1, .5]]), 2: numpy.array([[7, 1, .5]])}
  rngs = get_rngs(1)

  assert [{1: [Variation(POS=1, stop=2, REF='C', ALT='CAA', het=HOMOZYGOUS)],
           2: [Variation(POS=7, stop=8, REF='G', ALT='T', het=HET1),
               Variation(POS=17, stop=18, REF='G', ALT='T', het=HET2)]},
          {1: [Variation(POS=1, stop=2, REF='C', ALT='CAA', het=HET1)],
           2: [Variation(POS=7, stop=8, REF='G', ALT='T', het=HET2),
               Variation(POS=17, stop=18, REF='G', ALT='T', het=HET1)]}] == spawn(g1, g2, hot_spots, rngs)


def spawn_test2():
  """Spawn from parents (With denovo)."""
  ref = {
    1: 'CCTGACTGACTGACGTACGT',
    2: 'CCTGACGGACTGACGTGCGT'
  }

  g1 = {1: [Variation(1, 2, 'C', 'CAA', HET1)],
        2: [Variation(7, 8, 'G', 'T', HET2),
            Variation(17, 18, 'G', 'T', HET2)]}
  g2 = {1: [Variation(1, 2, 'C', 'CAA', HET2)],
        2: [Variation(7, 8, 'G', 'T', HET1),
            Variation(17, 18, 'G', 'T', HET1)]}

  # These hotspots ensure crossing over at these locii
  hot_spots = {1: numpy.array([[1, 1, .5]]), 2: numpy.array([[7, 1, .5]])}
  rngs = get_rngs(1)

  params_json = {
    "ss_variant_models": [
      {
        "snp": {
           "chromosome": [1],
           "phet": 1.0,
           "p": 0.1,
           "master_seed": 3
        }
      }
    ]
  }
  models = denovo.load_variant_models(params_json['ss_variant_models'])
  for m in models: reset_model(m)
  ch = spawn(g1, g2, hot_spots, rngs, 2, ref, models)
  correct_ch = [{1: [Variation(POS=1,stop=2,REF='C',ALT='CAA',het=HOMOZYGOUS),
                     Variation(POS=7,stop=8,REF='T',ALT='A',het=HET1),
                     Variation(POS=9,stop=10,REF='A',ALT='G',het=HET2)],
                 2: [Variation(POS=7,stop=8,REF='G',ALT='T',het=HET1),
                     Variation(POS=17,stop=18,REF='G',ALT='T',het=HET2)]},
                {1: [Variation(POS=1,stop=2,REF='C',ALT='CAA',het=HET1),
                     Variation(POS=7,stop=8,REF='T',ALT='G',het=HET1),
                     Variation(POS=16,stop=17,REF='T',ALT='A',het=HOMOZYGOUS)],
                 2: [Variation(POS=7,stop=8,REF='G',ALT='T',het=HET2),
                     Variation(POS=17,stop=18,REF='G',ALT='T',het=HET1)]}]
  assert correct_ch == ch, ch


def de_novo_population_test():
  """Spawn from parents (with denovo)"""
  params_json = {
    "denovo_variant_models": [
      {
        "snp": {
           "chromosome": [1],
           "phet": 0.9,
           "p": 0.1,
           "master_seed": 1
        }
      },
      {
        "snp": {
           "chromosome": [2],
           "phet": 0.9,
           "p": 0.1,
           "master_seed": 2
        }
      }
    ],
    "ss_variant_models": [
      {
        "snp": {
           "chromosome": [1],
           "phet": 1.0,
           "p": 0.075,
           "master_seed": 3
        }
      }
    ]
  }
  models = denovo.load_variant_models(params_json['denovo_variant_models'])
  for m in models: reset_model(m)
  ss_models = denovo.load_variant_models(params_json['ss_variant_models'])
  for m in ss_models: reset_model(m)

  ref = {
    1: 'CCTGACTGACTGACGTACGT',
    2: 'CCTGACGGACTGACGTGCGT'
  }

  pop = de_novo_population(ref, models, size=2)
  assert pop[0][1][0].POS == 16, pop
  assert pop[1][2][1].het == HOMOZYGOUS


@raises(AssertionError)
def one_generation_assert_test():
  """Too small a population: should not work"""
  one_generation([1], hot_spots={}, rngs={}, num_children_per_couple=2, ref=None, models=[])


def one_generation_test():
  """Create de novo population, then one set of children."""
  params_json = {
    "denovo_variant_models": [
      {
        "snp": {
           "chromosome": [1],
           "phet": 0.9,
           "p": 0.1,
           "master_seed": 1
        }
      },
      {
        "snp": {
           "chromosome": [2],
           "phet": 0.9,
           "p": 0.1,
           "master_seed": 2
        }
      }
    ],
    "ss_variant_models": [
      {
        "snp": {
           "chromosome": [1],
           "phet": 1.0,
           "p": 0.075,
           "master_seed": 3
        }
      }
    ]
  }
  models = denovo.load_variant_models(params_json['denovo_variant_models'])
  for m in models: reset_model(m)
  ss_models = denovo.load_variant_models(params_json['ss_variant_models'])
  for m in ss_models: reset_model(m)

  ref = {
    1: 'CCTGACTGACTGACGTACGT',
    2: 'CCTGACGGACTGACGTGCGT'
  }

  hot_spots = {1: numpy.array([[16, 1, .5]]), 2: numpy.array([[17, 1, 5]])}
  rngs = get_rngs(1)

  pop = de_novo_population(ref, models, size=2)
  children, parents = one_generation(pop, hot_spots=hot_spots, rngs=rngs, num_children_per_couple=2, ref=ref, models=ss_models)

  assert children[0][1][0].POS == 7, children
  assert children[1][2][0].REF == 'G'
  assert parents == [(1, 0)], parents


def population_simulation_test():
  """Create de novo population, then ten generations of descendants."""
  params_json = {
    "denovo_variant_models": [
      {
        "snp": {
           "chromosome": [1],
           "phet": 0.9,
           "p": 0.1,
           "master_seed": 1
        }
      },
      {
        "snp": {
           "chromosome": [2],
           "phet": 0.9,
           "p": 0.1,
           "master_seed": 2
        }
      }
    ],
    "ss_variant_models": [
      {
        "snp": {
           "chromosome": [1],
           "phet": 1.0,
           "p": 0.075,
           "master_seed": 3
        }
      }
    ]
  }
  models = denovo.load_variant_models(params_json['denovo_variant_models'])
  for m in models: reset_model(m)
  ss_models = denovo.load_variant_models(params_json['ss_variant_models'])
  for m in ss_models: reset_model(m)

  ref = {
    1: 'CCTGACTGACTGACGTACGT',
    2: 'CCTGACGGACTGACGTGCGT'
  }

  hot_spots = {1: numpy.array([[16, 1, .5]]), 2: numpy.array([[17, 1, 5]])}
  rngs = get_rngs(1)

  generations, parent_list = \
    population_simulation(ref, denovo_models=models, initial_size=2,
                          hot_spots=hot_spots, rngs=rngs, num_children_per_couple=2, ss_models=ss_models,
                          num_generations=10,
                          store_all_generations=True)

  assert generations[1][0][1][0].POS == 7
  assert generations[1][1][2][0].REF == 'G'
  assert parent_list[0] == [(1, 0)]
  assert len(generations) == 11