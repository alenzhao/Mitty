from mitty.lib.genome import FastaGenome
from mitty.population import founder_population
import mitty.denovo

from mitty.tests import example_fasta_genome


def founder_population_test():
  """Generating founder population with no denovo."""
  ancestor_models_json = [
    {
      "snp": {
         "chromosome": [1],
         "phet": 0.0,  # This is ignored for ancestors
         "p": 0.01,
      }
    },
    {
      "snp": {
         "chromosome": [2],
         "phet": 0.0,  # ignored for ancestors
         "p": 0.1,
      }
    }
  ]

  denovo_models_json = [
    {
      "snp": {
         "chromosome": [1],
         "phet": 0.0,  # This is ignored for ancestors
         "p": 0.0,
      }
    },
    {
      "snp": {
         "chromosome": [2],
         "phet": 0.0,  # ignored for ancestors
         "p": 0.0,
      }
    }
  ]

  ancestral_models = mitty.denovo.load_variant_model_list(ancestor_models_json)
  denovo_models = mitty.denovo.load_variant_model_list(denovo_models_json)
  ref = FastaGenome(seq_dir=example_fasta_genome, persist=True)

  pop, pool = founder_population(ancestral_models=ancestral_models,
                                 denovo_models=denovo_models,
                                 p_a=.95, p_het=0.5, ref=ref, pop_size=10, master_seed=1)

  for n in range(10):
    assert set(pop[n][1]) - set(pool[1]) == set()
    assert set(pop[n][2]) - set(pool[2]) == set()


# @nottest
# def chrom_crossover_test():
#   c1 = [
#       Variation(1, 2, 'C', 'CAA', HET_01),
#       Variation(3, 6, 'CAG', 'C', HET_10),
#       Variation(7, 8, 'G', 'T', HET_01),
#       Variation(9, 12, 'GTT', '', HET_10),
#       Variation(13, 16, 'GTT', 'TTG', HOMOZYGOUS),
#       Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS)
#   ]
#   crossover_idx = [1, 1, 1, 0, 1, 0]
#   correct_chrom = [
#       Variation(1, 2, 'C', 'CAA', HET_10),
#       Variation(3, 6, 'CAG', 'C', HET_01),
#       Variation(7, 8, 'G', 'T', HET_10),
#       Variation(9, 12, 'GTT', '', HET_10),
#       Variation(13, 16, 'GTT', 'TTG', HOMOZYGOUS),
#       Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS)
#   ]
#   c1_c = chrom_crossover(c1, crossover_idx)
#   assert correct_chrom == c1_c, c1_c
#   assert c1[0] == Variation(1, 2, 'C', 'CAA', HET_01)  # We shouldn't be changing our original
#
#
# @nottest
# def crossover_event_test():
#   c1 = [
#     Variation(1, 2, 'C', 'CAA', HET_01),
#     Variation(3, 6, 'CAG', 'C', HET_10),
#     Variation(7, 8, 'G', 'T', HET_01)
#   ]
#   c2 = [
#     Variation(2, 5, 'CAG', 'C', HET_10),
#     Variation(7, 8, 'G', 'T', HET_01)
#   ]
#   g1 = {1: c1, 2: c2}
#   crossover_idx = {1: [1, 0, 1], 2: [0, 1]}
#
#   correct_g = {
#     1: [
#       Variation(1, 2, 'C', 'CAA', HET_10),
#       Variation(3, 6, 'CAG', 'C', HET_10),
#       Variation(7, 8, 'G', 'T', HET_10)
#     ],
#     2: [
#       Variation(2, 5, 'CAG', 'C', HET_10),
#       Variation(7, 8, 'G', 'T', HET_10)
#     ]
#   }
#   assert correct_g == crossover_event(g1, crossover_idx)
#
#
# @nottest
# def pair_one_chrom_test():
#   """Merging: test c1 longer than c2"""
#   c1 = [
#     Variation(1, 2, 'C', 'CAA', HET_01),
#     Variation(3, 6, 'CAG', 'C', HET_10),
#     Variation(17, 18, 'G', 'T', HET_01),
#     Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS),
#   ]
#   c2 = [
#     Variation(7, 8, 'G', 'T', HET_01),
#   ]
#   which_copy = (1, 0)
#   correct_pairing = [
#     Variation(1, 2, 'C', 'CAA', HET_10),
#     Variation(17, 18, 'G', 'T', HET_10),
#     Variation(23, 26, 'GTT', 'TTG', HET_10),
#   ]
#   assert correct_pairing == pair_one_chrom(c1, c2, which_copy), pair_one_chrom(c1, c2, which_copy)
#
#
# @nottest
# def pair_one_chrom_test2():
#   """Merging: test c2 longer than c1"""
#   c1 = [
#     Variation(1, 2, 'C', 'CAA', HET_01)
#   ]
#   c2 = [
#     Variation(7, 8, 'G', 'T', HET_10),
#     Variation(17, 18, 'G', 'T', HET_10),
#   ]
#   which_copy = (1, 0)
#   correct_pairing = [
#     Variation(1, 2, 'C', 'CAA', HET_10),
#     Variation(7, 8, 'G', 'T', HET_01),
#     Variation(17, 18, 'G', 'T', HET_01),
#   ]
#   assert correct_pairing == pair_one_chrom(c1, c2, which_copy), pair_one_chrom(c1, c2, which_copy)
#
#
# @nottest
# def pair_one_chrom_test3():
#   """Merging: Comprehensive test"""
#   c1 = [
#     Variation(1, 2, 'C', 'CAA', HET_01),  # Tests homozygosity
#     Variation(3, 6, 'CAG', 'C', HET_10),  # Tests both variants are not on the copies chosen
#     Variation(17, 18, 'G', 'T', HET_01),  # Test zipper (several variants should come from c2 before we get to this)
#     Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS),  # Tests handling of homozygous variants
#     Variation(29, 30, 'T', 'G', HOMOZYGOUS)  # Tests unequal var list lengths
#   ]
#   c2 = [
#     Variation(1, 2, 'C', 'CAA', HET_10),
#     Variation(3, 6, 'CAG', 'C', HET_01),
#     Variation(7, 8, 'G', 'T', HET_01),
#     Variation(9, 12, 'GTT', '', HET_10),
#     Variation(15, 18, 'GTT', 'G', HET_10),  # Test partial overlap on different copies (should not affect each other)
#     Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS)
#   ]
#   which_copy = (1, 0)
#   correct_pairing = [
#     Variation(1, 2, 'C', 'CAA', HOMOZYGOUS),  # Tests homozygosity
#     Variation(9, 12, 'GTT', '', HET_01),
#     Variation(15, 18, 'GTT', 'G', HET_01),  # Test partial overlap on different copies (should not affect each other)
#     Variation(17, 18, 'G', 'T', HET_10),  # Test zipper (several variants should come from c2 before we get to this)
#     Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS),  # Tests handling of homozygous variants
#     Variation(29, 30, 'T', 'G', HET_10)  # Tests unequal var list lengths
#   ]
#   assert correct_pairing == pair_one_chrom(c1, c2, which_copy)
#
#
# @nottest
# def fertilize_one_test():
#   c11 = [
#     Variation(1, 2, 'C', 'CAA', HET_10)
#   ]
#   c12 = [
#     Variation(7, 8, 'G', 'T', HET_01),
#     Variation(17, 18, 'G', 'T', HET_01),
#   ]
#   g1 = {1: c11, 2: c12}
#
#   c21 = [
#     Variation(1, 2, 'C', 'CAA', HET_01)
#   ]
#   c22 = [
#     Variation(7, 8, 'G', 'T', HET_10),
#     Variation(17, 18, 'G', 'T', HET_10),
#   ]
#   g2 = {1: c21, 2: c22}
#   which_copy = {1: (0, 1), 2: (1, 0)}
#
#   correct_g3 = {
#     1: [
#       Variation(1, 2, 'C', 'CAA', HOMOZYGOUS)
#     ],
#     2: [
#       Variation(7, 8, 'G', 'T', HOMOZYGOUS),
#       Variation(17, 18, 'G', 'T', HOMOZYGOUS),
#     ]
#   }
#   assert correct_g3 == fertilize_one(g1, g2, which_copy)
#
#
# @nottest
# def place_crossovers_on_chrom_test():
#   """Cross over location generator"""
#   c1 = [
#     Variation(1, 2, 'C', 'CAA', HET_01),
#     Variation(3, 6, 'CAG', 'C', HET_10),
#     Variation(17, 18, 'G', 'T', HET_01),
#     Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS),
#     Variation(29, 30, 'T', 'G', HOMOZYGOUS)
#   ]
#
#   hot_spots = numpy.array([])  # No hotspots, no crossover
#   numpy.testing.assert_array_equal(numpy.array([0, 0, 0, 0, 0]),  # Hot spot is narrow and over first variant
#                                    place_crossovers_on_chrom(c1, hot_spots, numpy.random.RandomState(seed=1)))
#
#   hot_spots = numpy.array([[1, 1, .5]])
#   numpy.testing.assert_array_equal(numpy.array([1, 0, 0, 0, 0]),  # Hot spot is narrow and over first variant
#                                    place_crossovers_on_chrom(c1, hot_spots, numpy.random.RandomState(seed=1)))
#
#   hot_spots = numpy.array([[17, 1, .5]])
#   numpy.testing.assert_array_equal(numpy.array([0, 0, 1, 0, 0]),  # Hot spot is narrow and over third variant
#                                    place_crossovers_on_chrom(c1, hot_spots, numpy.random.RandomState(seed=1)))
#
#   hot_spots = numpy.array([[1, 1, .5], [17, 1, .5]])
#   numpy.testing.assert_array_equal(numpy.array([1, 0, 1, 0, 0]),  # Two narrow hot spots over first and third variants
#                                    place_crossovers_on_chrom(c1, hot_spots, numpy.random.RandomState(seed=1)))
#
#   hot_spots = numpy.array([[23, 1, 100], [17, 1, 100]])
#   numpy.testing.assert_array_equal(numpy.array([1, 1, 1, 1, 1]),  # Super broad hotspots, covers all
#                                    place_crossovers_on_chrom(c1, hot_spots, numpy.random.RandomState(seed=1)))
#
#   # A proper test is actually to do something like this (with c1 set as above)
#   # import pylab
#   # hot_spots = numpy.array([[3, 1, 2], [23, 1, 5]])
#   # rng = numpy.random.RandomState(seed=1)
#   # data = numpy.concatenate([numpy.nonzero(place_crossovers_on_chrom(c1, hot_spots, rng))[0] for n in range(100)])
#   # pylab.hist(data)
#   # pylab.plot()
#   # The denser your variant structure, the clearer is the sum of gaussian model
#
#
# @nottest
# def spawn_test():
#   """Spawn from parents (No denovo)."""
#   g1 = {1: [Variation(1, 2, 'C', 'CAA', HET_10)],
#         2: [Variation(7, 8, 'G', 'T', HET_01),
#             Variation(17, 18, 'G', 'T', HET_01)]}
#   g2 = {1: [Variation(1, 2, 'C', 'CAA', HET_01)],
#         2: [Variation(7, 8, 'G', 'T', HET_10),
#             Variation(17, 18, 'G', 'T', HET_10)]}
#
#   # These hotspots ensure crossing over at these locii
#   hot_spots = {1: numpy.array([[1, 1, .5]]), 2: numpy.array([[7, 1, .5]])}
#   rngs = get_rngs(1)
#
#   assert [{1: [Variation(POS=1, stop=2, REF='C', ALT='CAA', zygosity=HOMOZYGOUS)],
#            2: [Variation(POS=7, stop=8, REF='G', ALT='T', zygosity=HET_10),
#                Variation(POS=17, stop=18, REF='G', ALT='T', zygosity=HET_01)]},
#           {1: [Variation(POS=1, stop=2, REF='C', ALT='CAA', zygosity=HET_10)],
#            2: [Variation(POS=7, stop=8, REF='G', ALT='T', zygosity=HET_01),
#                Variation(POS=17, stop=18, REF='G', ALT='T', zygosity=HET_10)]}] == spawn(g1, g2, hot_spots, rngs)
#
#
# @nottest
# def spawn_test2():
#   """Spawn from parents (With denovo)."""
#   ref = {
#     1: 'CCTGACTGACTGACGTACGT',
#     2: 'CCTGACGGACTGACGTGCGT'
#   }
#
#   g1 = {1: [Variation(1, 2, 'C', 'CAA', HET_10)],
#         2: [Variation(7, 8, 'G', 'T', HET_01),
#             Variation(17, 18, 'G', 'T', HET_01)]}
#   g2 = {1: [Variation(1, 2, 'C', 'CAA', HET_01)],
#         2: [Variation(7, 8, 'G', 'T', HET_10),
#             Variation(17, 18, 'G', 'T', HET_10)]}
#
#   # These hotspots ensure crossing over at these locii
#   hot_spots = {1: numpy.array([[1, 1, .5]]), 2: numpy.array([[7, 1, .5]])}
#   rngs = get_rngs(1)
#
#   params_json = {
#     "ss_variant_models": [
#       {
#         "snp": {
#            "chromosome": [1],
#            "phet": 1.0,
#            "p": 0.1,
#            "master_seed": 3
#         }
#       }
#     ]
#   }
#   models = denovo.load_variant_models(params_json['ss_variant_models'])
#   for m in models: reset_model(m)
#   ch = spawn(g1, g2, hot_spots, rngs, 2, ref, models)
#   correct_ch = [{1: [Variation(POS=1,stop=2,REF='C',ALT='CAA',zygosity=HOMOZYGOUS),
#                      Variation(POS=7,stop=8,REF='T',ALT='A',zygosity=HET_10),
#                      Variation(POS=9,stop=10,REF='A',ALT='G',zygosity=HET_01)],
#                  2: [Variation(POS=7,stop=8,REF='G',ALT='T',zygosity=HET_10),
#                      Variation(POS=17,stop=18,REF='G',ALT='T',zygosity=HET_01)]},
#                 {1: [Variation(POS=1,stop=2,REF='C',ALT='CAA',zygosity=HET_10),
#                      Variation(POS=7,stop=8,REF='T',ALT='G',zygosity=HET_10),
#                      Variation(POS=16,stop=17,REF='T',ALT='A',zygosity=HOMOZYGOUS)],
#                  2: [Variation(POS=7,stop=8,REF='G',ALT='T',zygosity=HET_01),
#                      Variation(POS=17,stop=18,REF='G',ALT='T',zygosity=HET_10)]}]
#   assert_equal(correct_ch, ch)
#
#
# @nottest
# def de_novo_population_test():
#   """De novo population"""
#   params_json = {
#     "denovo_variant_models": [
#       {
#         "snp": {
#            "chromosome": [1],
#            "phet": 0.9,
#            "p": 0.1,
#            "master_seed": 1
#         }
#       },
#       {
#         "snp": {
#            "chromosome": [2],
#            "phet": 0.9,
#            "p": 0.1,
#            "master_seed": 2
#         }
#       }
#     ],
#     "ss_variant_models": [
#       {
#         "snp": {
#            "chromosome": [1],
#            "phet": 1.0,
#            "p": 0.075,
#            "master_seed": 3
#         }
#       }
#     ]
#   }
#   models = denovo.load_variant_models(params_json['denovo_variant_models'])
#   for m in models: reset_model(m)
#   ss_models = denovo.load_variant_models(params_json['ss_variant_models'])
#   for m in ss_models: reset_model(m)
#
#   ref = {
#     1: 'CCTGACTGACTGACGTACGT',
#     2: 'CCTGACGGACTGACGTGCGT'
#   }
#
#   pop = de_novo_population(ref, models, size=2)
#   assert pop[0][1][0].POS == 16, pop
#   assert pop[1][2][1].zygosity == HOMOZYGOUS
#
#
# @nottest
# @raises(AssertionError)
# def one_generation_assert_test():
#   """Too small a population: should not work"""
#   one_generation([1], hot_spots={}, rngs={}, num_children_per_couple=2, ref=None, models=[])
#
#
# @nottest
# def one_generation_test():
#   """Create de novo population, then one set of children."""
#   params_json = {
#     "denovo_variant_models": [
#       {
#         "snp": {
#            "chromosome": [1],
#            "phet": 0.9,
#            "p": 0.1,
#            "master_seed": 1
#         }
#       },
#       {
#         "snp": {
#            "chromosome": [2],
#            "phet": 0.9,
#            "p": 0.1,
#            "master_seed": 2
#         }
#       }
#     ],
#     "ss_variant_models": [
#       {
#         "snp": {
#            "chromosome": [1],
#            "phet": 1.0,
#            "p": 0.075,
#            "master_seed": 3
#         }
#       }
#     ]
#   }
#   models = denovo.load_variant_models(params_json['denovo_variant_models'])
#   for m in models: reset_model(m)
#   ss_models = denovo.load_variant_models(params_json['ss_variant_models'])
#   for m in ss_models: reset_model(m)
#
#   ref = {
#     1: 'CCTGACTGACTGACGTACGT',
#     2: 'CCTGACGGACTGACGTGCGT'
#   }
#
#   hot_spots = {1: numpy.array([[16, 1, .5]]), 2: numpy.array([[17, 1, 5]])}
#   rngs = get_rngs(1)
#
#   pop = de_novo_population(ref, models, size=2)
#   children, parents = one_generation(pop, hot_spots=hot_spots, rngs=rngs, num_children_per_couple=2, ref=ref, models=ss_models)
#
#   assert children[0][1][0].POS == 7, children
#   assert children[1][2][0].REF == 'G'
#   assert parents == [(1, 0)], parents
#
#
# @nottest
# def population_simulation_test():
#   """Create de novo population, then ten generations of descendants."""
#   params_json = {
#     "denovo_variant_models": [
#       {
#         "snp": {
#            "chromosome": [1],
#            "phet": 0.9,
#            "p": 0.1,
#            "master_seed": 1
#         }
#       },
#       {
#         "snp": {
#            "chromosome": [2],
#            "phet": 0.9,
#            "p": 0.1,
#            "master_seed": 2
#         }
#       }
#     ],
#     "ss_variant_models": [
#       {
#         "snp": {
#            "chromosome": [1],
#            "phet": 1.0,
#            "p": 0.075,
#            "master_seed": 3
#         }
#       }
#     ]
#   }
#   models = denovo.load_variant_models(params_json['denovo_variant_models'])
#   for m in models: reset_model(m)
#   ss_models = denovo.load_variant_models(params_json['ss_variant_models'])
#   for m in ss_models: reset_model(m)
#
#   ref = {
#     1: 'CCTGACTGACTGACGTACGT',
#     2: 'CCTGACGGACTGACGTGCGT'
#   }
#
#   hot_spots = {1: numpy.array([[16, 1, .5]]), 2: numpy.array([[17, 1, 5]])}
#   rngs = get_rngs(1)
#
#   import tempfile
#   pop_sim_data_dir = tempfile.mkdtemp()
#
#   parent_list, children = \
#     population_simulation(ref, denovo_models=models, initial_size=2,
#                           hot_spots=hot_spots, rngs=rngs, num_children_per_couple=2, ss_models=ss_models,
#                           num_generations=10,
#                           store_all_generations=True,
#                           vcf_prefix=os.path.join(pop_sim_data_dir, 'pop_sim_test'))
#
#   assert os.path.exists(os.path.join(pop_sim_data_dir, 'pop_sim_test_g0_p0.vcf'))
#   assert os.path.exists(os.path.join(pop_sim_data_dir, 'pop_sim_test_g1_p1.vcf.gz'))
#   assert os.path.exists(os.path.join(pop_sim_data_dir, 'pop_sim_test_g2_p0.vcf.gz.tbi'))
#   assert os.path.exists(os.path.join(pop_sim_data_dir, 'pop_sim_test_g9_p1.vcf'))
#   assert len(parent_list) == 10
#
#   import vcf
#   read_genome = parse_vcf(vcf.Reader(filename=os.path.join(pop_sim_data_dir, 'pop_sim_test_g10_p0.vcf.gz')), [1, 2])
#   assert children[0] == read_genome, [children[0], read_genome]
#
#   from shutil import rmtree
#   rmtree(pop_sim_data_dir)  # Clean up after ourselves