"""This is to simulate sexual reproduction in a population of individuals complete with recombination (crossover).

In contrast to denovo, generations does not directly write out VCF files for each sample. The variant data are stored in
a sqlite database


The general model is to, first, generate a founder population of completely unrelated samples. Pairs of samples are then
picked for mating and children are generated. If the number of children exceeds the maximum population size we prune
away the children based on fitness values. A genealogy tree of the pairings are kept. Population breaks, where
sub-populations are isolated at particular generations can also be modeled. This main overall behavior is shaped by
different plugin models that describe the details of each stage.


The stock plugins are described later below

1. mating plugin: randomly picks pairs of mates from the pool. A picked pair is rejected if they have the same parent
N generations back. The number of children is determined by the combined fitness values of the parents

2. cross-over: has a flat probability of cross over points. If you find hotspots, your algorithm is flawed

3. fitness: has an (arbitrary) map of regions where mutations are deleterious and where they are beneficial and where
mutations are recessive and dominant

4. fracture: break the given population into two separate parts


NOTES:

Should be able to start from the middle of a simulation - basically, we send the same file in and mention a generation
in the sim and we should be able to pickup from there, given the descendant tree, random number seed sequence and the
full genomes of the last generation are available. We should, similarly, be able to take an arbitrary, existing
population as the parent list and move the simulation forward.
"""
import random


class Sample:
  """This represents one genome sample.

  Attributes:
    generation - which generation is this sample
    serial     - what number in this generation
    name       - string for the form 'generation:serial'
    parents    - list of references to parent Samples
    genome     - the genome dictionary of this sample
    fitness    - computed fitness value of the genome
  """
  def __init__(self, generation, serial):
    self.name = '{:d}:{:d}'.format(generation, serial)  # Our unique sample name
    self.generation, self.serial, self.parents = generation, serial, []
    self.genome, self.fitness = {}, 0  # Fitness is filled in with the genome

  def __hash__(self):
    return self.name.__hash__()

  def __eq__(self, other):
    return self.__hash__() == other.__hash__()

  def __repr__(self):
    return '{:s} ({:d})'.format(self.name, self.fitness)


def incestuous_mating(s1, s2, incest_generations=2):
  """Return False if the samples are not related up to the generations indicated. (Find if there is an LCA within x levels)
  :param s1, s2: samples to test
  :param incest_generations: number of generations to go back for relatedness test
  :returns boolean
  """
  parent_list = [p for s in [s1, s2] for p in s.parents]
  for n in range(incest_generations):
    if len(set(parent_list)) < len(parent_list):
      return True
    parent_list = [p for s in parent_list for p in s.parents]
  return False


def mate_us(parent_list, incest_generations=1):
  """Simple mating generator pops pairs of parents from the list until the list is exhausted.
  :param parent_list: list of Samples that will be parents. Items will be popped from this list, just so you know
  :param incest_generations: How many generations back should we go to check for incest
  :returns mating_list - list of pairs of legal mates
           incest_list - two lists of samples. Corresponding samples in the lists are related
  """
  mating_list = []
  incest_list = [[], []]  # The folks we couldn't mate, kept in separate lists to increase their chances of mating
  while len(parent_list) > 1:
    s1 = parent_list.pop()
    s2 = parent_list.pop()
    if incestuous_mating(s1, s2, incest_generations):
      incest_list[0] += [s1]
      incest_list[1] += [s2]
      continue
    mating_list += [[s1, s2]]
  return mating_list, incest_list


def generate_matings(parent_list, incest_generations=1, master_seed=2):
  """Simple mating generator that picks pairs of parents randomly until the pool is exhausted.
  :param parent_list: list of Samples that will be parents. This will be shuffled in place and items will be popped
                      from it, just so you know
  :param incest_generations: How many generations back should we go to check for incest
  :param master_seed: Integer seed for RNG

  We create a list of

  """
  rnd = random.Random(master_seed)
  rnd.shuffle(parent_list)
  mating_list, incest_list = mate_us(parent_list, incest_generations)

  # Make an attempt on the rejects
  mating_list2, _ = mate_us(incest_list[0], incest_generations)
  mating_list3, _ = mate_us(incest_list[1], incest_generations)

  return mating_list + mating_list2 + mating_list3


def generate_pair():
  pass


#
#
# def chrom_crossover(c1, crossover_idx):
#   """cross_over_idx is a list the same size as c1 with a 1 (indicating a crossover should be done) or 0 (no crossover).
#   """
#   c2 = []
#   for c, idx in zip(c1, crossover_idx):
#     new_c = vcopy(c)  # Valid for no crossover or HOM
#     if idx == 1 and c.zygosity != HOM:
#       if c.zygosity == HET_10:
#         new_c.zygosity = HET_01
#       else:
#         new_c.zygosity = HET_10
#     c2 += [new_c]
#   return c2
#
#
# def crossover_event(g1, crossover_idx):
#   return {chrom: chrom_crossover(g1[chrom], crossover_idx[chrom]) for chrom in g1}
#
#
# def pair_one_chrom(c1, c2, which_copy):
#   """Has a resemblance to merge sort BUT WITH GENOMIC VARIANTS!
#   Parameters
#   ----------
#   c1, c2     : tuple list
#                (POS, POS2, REF, ALT, HET)
#   which_copy : tuple
#                e.g. (0, 1) telling us which chromosome of each individual to take
#
#   Returns
#   -------
#   c3         : tuple list
#                (POS, POS2, REF, ALT, HET)
#   """
#   c1_iter, c2_iter = c1.__iter__(), c2.__iter__()
#   c3 = []
#   # Try the zipper
#   l1, l2 = next(c1_iter, None), next(c2_iter, None)
#   while l1 is not None and l2 is not None:
#
#     if (l1.zygosity == HET_10 and which_copy[0] == 1) or (l1.zygosity == HET_01 and which_copy[0] == 0):
#       l1 = next(c1_iter, None)
#       continue
#
#     if (l2.zygosity == HET_10 and which_copy[1] == 1) or (l2.zygosity == HET_01 and which_copy[1] == 0):
#       l2 = next(c2_iter, None)
#       continue
#
#     if vcopy(l1, zygosity=HOM) == vcopy(l2, zygosity=HOM):  # Homozygous
#       c3 += [vcopy(l1, zygosity=HOM)]
#       l1, l2 = next(c1_iter, None), next(c2_iter, None)
#       continue
#
#     if l1.POS <= l2.POS:
#       c3 += [vcopy(l1, zygosity=HET_10)]
#       l1 = next(c1_iter, None)
#     else:
#       c3 += [vcopy(l2, zygosity=HET_01)]
#       l2 = next(c2_iter, None)
#
#   # Now pick up any slack
#   while l1 is not None:
#     if (l1.zygosity == HOM) or (l1.zygosity == HET_10 and which_copy[0] == 0) or (l1.zygosity == HET_01 and which_copy[0] == 1):
#       c3 += [vcopy(l1, zygosity=HET_10)]
#     l1 = next(c1_iter, None)
#
#   while l2 is not None:
#     if (l2.zygosity == HOM) or (l2.zygosity == HET_10 and which_copy[1] == 0) or (l2.zygosity == HET_01 and which_copy[1] == 1):
#       c3 += [vcopy(l2, zygosity=HET_01)]
#     l2 = next(c2_iter, None)
#
#   return c3
#
#
# def fertilize_one(g1, g2, which_copy):
#   #return {chrom: pair_one_chrom(c1, c2, wc) for (chrom, c1), (_, c2), wc in zip(g1.iteritems(), g2.iteritems(), which_copy)}
#   return {chrom: pair_one_chrom(g1[chrom], g2[chrom], which_copy[chrom]) for chrom in g1}
#
#
# def place_crossovers_on_chrom(c1, hot_spots, rng):
#   """Stock cross over generator. Place segments based on gaussian distribution around hot_spots
#   The model is a mixture-gaussian model: the hotspot has an amplitude and a width (sd of the gaussian). The amplitude of
#   the gaussian at a certain locus determines the probability of the variants there crossing over. The probability of
#   crossover at any point is the sum of all Gaussians clipped at +1.0
#
#   Parameters
#   ----------
#   c1               : list of tuples
#                      [(start, stop, <any other payload>) ...]. Should be sorted
#   hot_spots        : numpy array
#                      n x 3
#                      [(center, max_p, width) ... ] 0 <= max_p <= 1
#   rng              : random number generator
#
#   Returns
#   -------
#   idx              : numpy array
#                      As used by chrom_crossover
#   """
#   if hot_spots.shape[0] == 0:
#     return [0] * len(c1)  # No hotspots, no crossovers
#   #x = numpy.array(c1, dtype=[('st', float), ('a', 'c'), ('b', 'c'), ('c', 'c'), ('d', 'c')])['st']
#   x = [v.POS for v in c1]
#   X, C = numpy.meshgrid(x, hot_spots[:, 0])
#   _, A = numpy.meshgrid(x, hot_spots[:, 1])
#   _, W = numpy.meshgrid(x, hot_spots[:, 2])
#   p = (A * numpy.exp(-(((X - C) / W) ** 2))).sum(axis=0).clip(0, 1.0)
#   #return numpy.nonzero(rng.uniform(size=p.size) < p)[0]
#   return rng.uniform(size=p.size) < p
#
#
# def place_crossovers(g1, hot_spots, rng):
#   return {chrom: place_crossovers_on_chrom(g1[chrom], hot_spots[chrom], rng) for chrom in g1}
#
#
# def which_copies(g1, rng):
#   """Simple 50/50 coin toss to decide which chromosome copy gets picked during fertilization."""
#   return {chrom: wc for chrom, wc in zip(g1, rng.randint(0, 2, size=(len(g1), 2)))}
#
#
# def get_rngs(seed):
#   """Creates independent RNGs for our random processes from the master seed given."""
#   rng_names = ['cross_over', 'chrom_copy', 'couple_chose']
#   return {k: numpy.random.RandomState(seed=sub_seed)
#           for k,sub_seed in zip(rng_names, numpy.random.RandomState(seed=seed).randint(100000000, size=len(rng_names)))}
#
#
# def spawn(g1, g2, hot_spots={}, rngs={}, num_children=2, ref=None, models=[]):
#   """If you don't want denovo mutations don't send in ref_fp and models"""
#   if g1.keys() != g2.keys():
#     raise RuntimeError('Two genomes have unequal chromosomes')
#   children = []
#   for _ in range(num_children):
#     g1_cross = crossover_event(g1, place_crossovers(g1, hot_spots, rngs['cross_over']))
#     g2_cross = crossover_event(g2, place_crossovers(g2, hot_spots, rngs['cross_over']))
#     if ref is not None and len(models) > 0:  # We want denovo mutations
#       g1_cross = mitty.denovo.apply_variant_models_to_genome(g1_cross, ref, models)
#       g2_cross = mitty.denovo.apply_variant_models_to_genome(g2_cross, ref, models)
#     children.append(fertilize_one(g1_cross, g2_cross, which_copies(g1, rngs['chrom_copy'])))
#   return children
#
#
# # def find_mates(pop, related_level=0, parent_list=[]):
# #   """Related level is the number of generations we need to go back to ensure they are not related."""
# #   our_pop = list(pop)
# #   pairs = []
# #   n = 0
# #   while n < len(our_pop):
# #     if our_pop[n] is None:
# #       n += 1
# #       continue
# #     m = 0
# #     while m < len(pop):
# #
# #     this_pair = [n]
#
#
# def one_generation(pop, hot_spots={}, rngs={}, num_children_per_couple=2, ref=None, models=[]):
#   """We take pairs from the population without replacement, cross over the chromosomes, sprinkle in denovo mutations
#   then shuffle the chromosomes during fertilization. If the population size is odd, we discard one parent"""
#   assert len(pop) >= 2, 'Need at least two parents'
#   mates = rngs['couple_chose'].permutation(len(pop))
#   children = []
#   parents = []
#   for n in range(2*(len(mates)/2))[::2]:  #Todo fix handling of odd sizes
#     children += spawn(pop[mates[n]], pop[mates[n + 1]], hot_spots=hot_spots, rngs=rngs, num_children=num_children_per_couple,
#                       ref=ref, models=models)
#     parents += [(mates[n], mates[n + 1])]
#   return children, parents
#
#
# def population_simulation(ref, denovo_models=[], initial_size=10,
#                           hot_spots={}, rngs={}, num_children_per_couple=2, ss_models=[],
#                           num_generations=10,
#                           store_all_generations=True,
#                           vcf_prefix=''):
#   # The initial population
#   pop = de_novo_population(ref, denovo_models, size=initial_size)
#   logger.debug('Created de novo population of {:d} individuals'.format(len(pop)))
#   save_one_generation(pop, generation=0, vcf_prefix=vcf_prefix)
#
#   parent_list = []
#   # The generational loop
#   this_gen = children = pop
#   for n in range(1, num_generations + 1):
#     children, parents = one_generation(this_gen, hot_spots=hot_spots, rngs=rngs,
#                                        num_children_per_couple=num_children_per_couple,
#                                        ref=ref, models=ss_models)
#     logger.debug('Created generation {:d} ({:d})'.format(n, len(children)))
#     parent_list.append(parents)
#     if store_all_generations:
#       save_one_generation(children, generation=n, vcf_prefix=vcf_prefix)
#     this_gen = children
#
#   if not store_all_generations:
#       save_one_generation(pop, generation=n, vcf_prefix=vcf_prefix)
#
#   return parent_list, children
#
#
# def save_one_generation(pop, generation=0, vcf_prefix='pop'):
#   """Given a population of genomes same them as a collection of vcf files with a prefix and sequential numbering."""
#   for n, p in enumerate(pop):
#     vcf_save_gz(p, vcf_gz_name='{:s}_g{:d}_p{:d}.vcf.gz'.format(vcf_prefix, generation, n),
#                 sample_name='g{:d}p{:d}'.format(generation, n))
#
#
# def save_multiple_generations(gen_pop, vcf_prefix='pop'):
#   for gen, pop in enumerate(gen_pop):
#     save_one_generation(pop, generation=gen, vcf_prefix=vcf_prefix)
#
#
# def main(ref_fname='../../Data/ashbya_gossypii/ashbya.h5', hot_spots={}, params_json={}, master_seed=None,
#          initial_size=10, num_children_per_couple=2, num_generations=20,
#          store_all_generations=True,
#          out_dir='out',
#          out_prefix='pop'):
#   import mitty.denovo as denovo
#   import shutil
#   import os
#   shutil.rmtree(out_dir)
#   os.makedirs(out_dir)
#
#   models = denovo.load_variant_models(params_json['denovo_variant_models'])
#   if master_seed is not None:
#     denovo.apply_master_seed(models, master_seed=int(master_seed))
#     rngs = get_rngs(seed=int(master_seed))
#   else:
#     rngs = get_rngs(seed=1)
#   ss_models = denovo.load_variant_models(params_json['ss_variant_models'])
#   ref = genome.FastaGenome(ref_fname)
#   parent_list, _ = \
#     population_simulation(ref, denovo_models=models, initial_size=initial_size,
#                           hot_spots=hot_spots, rngs=rngs, num_children_per_couple=num_children_per_couple,
#                           ss_models=ss_models,
#                           num_generations=num_generations,
#                           store_all_generations=store_all_generations,
#                           vcf_prefix=os.path.join(out_dir, out_prefix))
#
# if __name__ == "__main__":
#   import json
#   import docopt
#
#   if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
#     docopt.docopt(__doc__, ['-h'])
#   else:
#     args = docopt.docopt(__doc__, version=__version__)
#   if args['plugins']:
#     mitty.denovo.print_plugin_list()
#     exit(0)
#   if args['explain']:
#     mitty.denovo.explain_plugin(args['<plugin>'])
#     exit(0)
#
#   logging.basicConfig(level=logging.WARNING)
#   if args['-v']:
#     logger.setLevel(logging.DEBUG)
#
#   hs = json.load(open(args['--hs'], 'r'))
#   hot_spots = {int(c): numpy.array(v) for c, v in hs.iteritems()}
#   pj = json.load(open(args['--paramfile'], 'r'))
#   main(ref_fname=args['--wg'], hot_spots=hot_spots, params_json=pj, master_seed=args['--master_seed'],
#        initial_size=int(args['--init_size']), num_children_per_couple=int(args['--children']),
#        num_generations=int(args['--gens']), store_all_generations=False,
#        out_dir=args['--outdir'], out_prefix=args['--outprefix'])

