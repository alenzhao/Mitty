"""This module defines a structure to carry variation information and some functions to perform operations on
genomes (collections of variations). For a description of design choices and algorithms please see the developer
documentation
"""

# Types of zygosity
ABSENT = 0
HOMOZYGOUS = 3
HET_01 = 1
HET_10 = 2
#       00     01     10     11
GT = ['0|0', '0|1', '1|0', '1|1']  # This needs to match rev 3 definitions


cdef class VariationData:
  """VD(pos, stop, REF, ALT)
  Attributes:
    POS
    footprint
    REF
    ALT"""
  cdef public:
    long int POS, stop
    str REF, ALT

  def __cinit__(self, int _pos=0, int _stop=0, str _ref='', str _alt=''):
    self.POS, self.stop, self.REF, self.ALT = _pos, _stop, _ref, _alt


def new_variation(int _pos=0, int _stop=0, str _ref='', str _alt='', char het=0, char recessive=0, float fitness=0):
  return Variation(VariationData(_pos, _stop, _ref, _alt), het, recessive, fitness)


cdef class Variation:
  """Carries *reference* to a VariationData instance and heterozygosity, recessiveness and fitness data."""
  cdef public:
    VariationData vd
    unsigned char het, recessive, fitness

  # http://docs.cython.org/src/userguide/special_methods.html
  def __cinit__(self, VariationData _vd, unsigned char het=HOMOZYGOUS, unsigned char rec=1, unsigned char fit=128):
    self.vd, self.het, self.recessive, self.fitness = _vd, het, rec, fit


  # Need this for set to work (https://docs.python.org/2/reference/datamodel.html#object.__hash__)
  def  __hash__(self):
    return self.vd.POS ^ self.vd.REF.__hash__() ^ self.vd.ALT.__hash__()


  cpdef bint eq(self, Variation other):
    # We don't consider zygosity, fitness or recessiveness
    return self.vd.POS == other.vd.POS and self.vd.stop == other.vd.stop and \
           self.vd.REF == other.vd.REF and self.vd.ALT == other.vd.ALT

  # http://docs.cython.org/src/userguide/special_methods.html
  def __richcmp__(self, Variation other, int op):
    if op == 2:
      return self.eq(other)
    elif op == 3:
      return not self.eq(other)
    else:
      return False

  def __repr__(self):
    """Return a human friendly printed representation of the variant."""
    return '(POS={0},stop={1},REF={2},ALT={3},zygosity={4},r={5},fit={6})'.\
      format(self.vd.POS, self.vd.stop, self.vd.REF, self.vd.ALT, GT[self.het], self.recessive, self.fitness)


def copy_variant_sequence(list c1, idx=None, het=None):
  """copy_variant_sequence(list c1, idx, zygosity)
  :param list c1: a sequence of variants that need to be copied over.
  :rtype generator: Generated copy of the variant list. Note that we always share the VariantData"""
  cdef Variation v1
  cdef int n
  if idx is None:
    return (Variation(v1.vd, v1.het, v1.recessive, v1.fitness) for v1 in c1)
  elif het is None:
    return (Variation(c1[n].vd, c1[n].het, c1[n].recessive, c1[n].fitness) for n in idx)
  else:
    return (Variation(c1[n].vd, het[n2], c1[n].recessive, c1[n].fitness) for n2, n in enumerate(idx))


cdef inline Variation copy_variant(Variation v1):
  return Variation(v1.vd, v1.het, v1.recessive, v1.fitness)


def merge_variants(c1, c2):
  """merge_variants(c1, c2)
  Given an existing chromosome (list of variants) merge a new list of variants into it in zipper fashion. c1 has
  priority (collisions are resolved in favor of c1)

  :param c1: The original variants (iterable)
  :param c2: The proposed new variants (iterable)
  :rtype list: The resultant variant list with collisions arbitrated

  Algorithm::

    o(x,y) = True if x and y overlap
           = False otherwise
    c1 = existing list of variants
    c2 = denovo list of variants
    c3 = new list of variants being built

    o(c1, c2) = True: add(c1) c1++, c2++
              = False:
                c1 < c2 ? add(c1) c1++
                else:
                  o(c3, c2) = True: c2++
                            = False: add(c2), c2++

  """
  return c_merge_variants(c1, c2)


cdef inline bint overlap(Variation x, Variation y):
  # Returns true if the footprints of the variations overlap.
  if x is None or y is None: return False
  if y.vd.POS - 1 <= x.vd.POS <= y.vd.stop + 1 or y.vd.POS - 1 <= x.vd.stop <= y.vd.stop + 1 or \
                  x.vd.POS <= y.vd.POS - 1 <= y.vd.stop + 1 <= x.vd.stop:  # Potential overlap
    if x.het == y.het or x.het == HOMOZYGOUS or y.het == HOMOZYGOUS:  # Definite overlap
      return True
  return False


#TODO: refactor this to be faster? Not use copy_variant? Code can be made cleaner
cdef c_merge_variants(list c1, list c2):
  cdef:
    int n1_max = len(c1), n2_max = len(c2)
    int n1 = 0, n2 = 0, n3 = 0
  c3 = [None] * (n1_max + n2_max)  # This is the maximum size of c3
  while n1 < n1_max and n2 < n2_max:
    if overlap(c1[n1], c2[n2]): # This will collide. Advance c2 and redo. Fixes case for merge_test8
      n2 += 1
    else:  # No collision
      if c1[n1].vd.POS <= c2[n2].vd.POS:  # Zip-in c1 as it comes before c2
        c3[n3] = copy_variant(c1[n1])
        n1 += 1
        n3 += 1
      else:  # Zip-in c2 as it comes before next c1
        if n3==0 or not overlap(c3[n3 - 1], c2[n2]):
          c3[n3] = copy_variant(c2[n2])
          n3 += 1
        n2 +=1  # Need to advance c2 anyway

  # Now copy over slack
  while n1 < n1_max:  # Smooth sailing, just copy over the rest of the original (c1)
    c3[n3] = copy_variant(c1[n1])
    n1 += 1
    n3 += 1

  while n2 < n2_max:  # Need to test each new (c2) for clashes with itself
    if n3==0 or not overlap(c3[n3 - 1], c2[n2]):
      c3[n3] = copy_variant(c2[n2])
      n3 += 1
    n2 += 1  # Need to advance c2 anyway

  return c3[:n3]


def copy_missing_chromosomes(g1, g2):
  """Copy any chromosomes found in g2 but not in g1 onto g1. g1 is changed in place

  :param dict g1: genome
  :param dict g2: genome
  :returns: changes g1 in place"""
  missing_chrom = set(g2.keys()) - set(g1.keys())
  cdef int ch
  for ch in missing_chrom:
    g1[ch] = list(copy_variant_sequence(g2[ch]))


def merge_genomes(g1={}, g2={}):
  """Given two genomes run merge_variants(c1, c2) on each chromosome. g1 has priority. In the end copy any chromosomes
  present in g2 but not in g1 into g1 (so that g1 is a superset of g2)

  :param dict g1: genome
  :param dict g2: genome
  :returns: new genome
  """
  g3 = {}
  for chrom, c2 in g2.iteritems():
    g3[chrom] = merge_variants(g1.get(chrom, []), c2)
  copy_missing_chromosomes(g3, g1)  # The previous loop merges all chr in c2.
                                    # Now, we need to consider all chr in c1 but not in c2
  return g3
