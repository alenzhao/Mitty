"""This module defines a structure to carry variation information and some functions to perform operations on
genomes (collections of variations). For a description of design choices and algorithms please see the developer
documentation
"""
import itertools  # Because zip does not play well with iterators

# Types of zygosity
cpdef enum:
  ABSENT = 0
  HET_01 = 1
  HET_10 = 2
  HOM = 3
#       00     01     10     11
GT = ['0|0', '0|1', '1|0', '1|1']  # This needs to match rev 3 definitions


cdef class Variant:
  """Variant(pos, stop, REF, ALT)
  A lightweight class to carry a single data

  Attributes:
    pos    - position of data
    stop   - where does footprint of data on reference (position of last base + 1 in the REF entry)
    hash   - used to check equality between variants
    ref    - reference sequence
    alt    - alternate sequence

  Each Variant is placed in a master list dictionary. The index in the dictionary is computed from the pos and the
  number allele this is at that position. This index is then used as the rowid when we save the data to the database.
  """
  cdef public:
    unsigned long pos, stop, hash, index
    bytes ref, alt

  def __cinit__(self, unsigned long pos, unsigned long stop, bytes ref, bytes alt):
    self.pos, self.stop, self.ref, self.alt = pos, stop, ref, alt
    self.hash = self.pos ^ hash(self.ref) ^ hash(self.alt)

  def __repr__(self):
    return '({:d}, {:d}, {:s}, {:s})'.format(self.pos, self.stop, self.ref, self.alt)

  def as_tuple(self):
    return self.pos, self.stop, self.ref, self.alt


cpdef Variant add_novel_variant_to_master(Variant v, dict ml):
  """Add this new variant to the dictionary (master list). Intelligently handle case where we have an existing allele
  at the same locus.
  :param v: This variant in the form of a Variation instance
  :param ml: Python dictionary storing all variants in a chromosome
  :returns v or existing variant in master list to which this variant is identical
  """
  cdef unsigned short n = 0
  cdef unsigned long index = v.pos << 16 | n
  while index in ml:  # Find us an empty spot
    if ml[index].hash == v.hash:
      return ml[index] # This is identical to an existing variant
    n += 1
    index = v.pos << 16 | n
  ml[index] = v
  v.index = index
  return v


cdef class SampleVariant:
  """Represents a data in a sample. It points to the Variant information and carries genotype info. It also has a
  pointer to the next SampleVariant so we can chain it together into a Sample

  SampleVariant(gt, Variant)
  Attributes:
    gt      - zygosity (genotype) information
    data - (pointer to) the Variant in the master list
    """
  cdef public:
    unsigned char gt
    Variant data
    SampleVariant next

  def __cinit__(self, unsigned char gt, Variant variant):
    self.gt, self.data, self.next = gt, variant, None

  def __repr__(self):
    return '{:d} {:s}'.format(self.data.pos, GT[self.gt])


def create_sample_iterable(pos, ref, alt, gt):
  """Given lists/generators of the bare variant data (e.g. from a plugin) return us an iterator that will produce a
  stream of individual SampleVariants."""
  for p, r, a, g in itertools.izip(pos, ref, alt, gt):
    yield SampleVariant(g, Variant(p, p + len(r), r, a))


cdef class Sample:
  """Represents a sample as a linked list. We always go sequentially through a sample, and for denovo variants
  we may add samples into the middle of the sample which is a good fit for a ll

  Rules:
    There is a dummy node, which head (and, initially, tail) point to.
    appends attach to tail.next and reposition tail. This is how we create the list
    inserts attach in between cursor cursor.next
    advance advances cursor and returns cursor.next
  """
  cdef public:
    str label
    SampleVariant head, cursor, tail
    unsigned long length

  def __cinit__(self, str label):
    self.head = self.tail = SampleVariant(ABSENT, None)  # Dummy node
    self.label, self.cursor, self.length = label, None, 0

  def __len__(self):
    return self.length

  def __iter__(self):
    return SampleIterator(self)

  def to_list(self):
    return [sv for sv in self]

  cpdef rewind_cursor(self):
    self.cursor = None

  cpdef SampleVariant advance(self):
    if self.cursor is None:  # We are at the head of the list
      self.cursor = self.head
    else:
      if self.cursor.next is not None:  # OK to advance
        self.cursor = self.cursor.next
    return self.cursor.next  # Will be None for empty list

  cpdef insert(self, SampleVariant sv):
    """Place sv right after cursor. cursor will point to original cursor.next after this. In practice, this places
    sv right before the variant advance just spit out"""
    sv.next = self.cursor.next
    self.cursor.next = sv
    self.cursor = sv
    if sv.next is None:
      self.tail = sv
    self.length += 1

  cpdef append(self, SampleVariant sv):
    """This is how we grow the list."""
    self.tail.next = sv
    self.tail = sv
    self.length += 1


cdef class SampleIterator:
  """A class that lets us iterate over the sample."""
  cdef SampleVariant this

  def __cinit__(self, Sample s):
    self.this = s.head.next  # Head is a dummy!

  def __next__(self):
    if self.this is None:
      raise StopIteration()
      return None
    else:
      result = self.this
      self.this = self.this.next
      return result


cdef inline bint overlap(SampleVariant v1, SampleVariant v2):
  """Returns true if the footprints of the variations overlap. This is used when applying denovo mutations to a genome
  :param v1: Variant from original sample
  :param v2: Proposed Variant
  :param gt: proposed genotype for this sample
  :returns Bool
  """
  if v2.data.pos - 1 <= v1.data.pos <= v2.data.stop + 1 or v2.data.pos - 1 <= v1.data.stop <= v2.data.stop + 1 or \
     v1.data.pos <= v2.data.pos - 1 <= v2.data.stop + 1 <= v1.data.stop:  # Potential overlap
    if v1.gt == v2.gt or v1.gt == HOM or v2.gt == HOM:  # Definite overlap
      return True
  return False


cpdef add_denovo_variants_to_sample(Sample s, dnv, dict ml):
  """add_denovo_variants_to_sample(s, dnv, ml)
  Given an existing Sample s (list of variants) merge a new list of variants (dnv) into it in zipper fashion. s has
  priority (collisions are resolved in favor of s). As new variants are accepted, add them to the master list

  :param c1: The original variants
  :param c2: The proposed new variants (iterator, convenient to use create_sample_iterable)

  s and ml are modified in place
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
  cdef:
    SampleVariant s1 = s.advance(), s2 = next(dnv, None)
  while s1 is not None and s2 is not None:
    if overlap(s1, s2): # This will collide. Advance dnv and redo. Fixes case for merge_test8
      s2 = next(dnv, None)
    else:  # No collision
      if s1.data.pos <= s2.data.pos:  # Advance s until we come to s2
        s1 = s.advance()
      else:  #s1 is past s2 and there is no collision, good to add
        s2.data = add_novel_variant_to_master(s2.data, ml)  # If this denovo already exists in the master list, use that
        s.insert(s2)
        s2 = next(dnv, None)

  while s2 is not None:  # All the remaining denovo variants can just be appended to the sample
    s2.data = add_novel_variant_to_master(s2.data, ml)
    s.append(s2)
    s2 = next(dnv, None)


# cpdef merge_variants(list c1, list c2):
#   """merge_variants(c1, c2)
#   Given an existing chromosome (list of variants) merge a new list of variants into it in zipper fashion. c1 has
#   priority (collisions are resolved in favor of c1)
#
#   :param c1: The original variants (iterable)
#   :param c2: The proposed new variants (iterable)
#   :rtype list: The resultant data list with collisions arbitrated
#
#   Algorithm::
#
#     o(x,y) = True if x and y overlap
#            = False otherwise
#     c1 = existing list of variants
#     c2 = denovo list of variants
#     c3 = new list of variants being built
#
#     o(c1, c2) = True: add(c1) c1++, c2++
#               = False:
#                 c1 < c2 ? add(c1) c1++
#                 else:
#                   o(c3, c2) = True: c2++
#                             = False: add(c2), c2++
#
#   """
#   cdef:
#     int n1_max = len(c1), n2_max = len(c2)
#     int n1 = 0, n2 = 0, n3 = 0
#     Variant v1, v2
#     Genotype g1, g2
#     list c3 = [None] * (n1_max + n2_max)  # This is the maximum size of c3
#   while n1 < n1_max and n2 < n2_max:
#     g1, g2 = c1[n1], c2[n2]
#     v1, v2 = ml[g1.index], ml[g2.index]
#     if overlap(v1, v2): # This will collide. Advance c2 and redo. Fixes case for merge_test8
#       n2 += 1
#     else:  # No collision
#       if v1.pos <= v2.pos:  # Zip-in c1 as it comes before c2
#         c3[n3] = g1
#         n1 += 1
#         n3 += 1
#       else:  # Zip-in c2 as it comes before next c1
#         if n3==0 or not overlap(c3[n3 - 1], v2):
#           c3[n3] = v2
#           n3 += 1
#         n2 +=1  # Need to advance c2 anyway
#
#   # Now copy over slack
#   while n1 < n1_max:  # Smooth sailing, just copy over the rest of the original (c1)
#     c3[n3] = v1
#     n1 += 1
#     n3 += 1
#
#   while n2 < n2_max:  # Need to test each new (c2) for clashes with itself
#     if n3==0 or not overlap(c3[n3 - 1], v2):
#       c3[n3] = c2[n2]
#       n3 += 1
#     n2 += 1  # Need to advance c2 anyway
#
#   return c3[:n3]


# def copy_missing_chromosomes(dict g1, dict g2):
#   """Copy any chromosomes found in g2 but not in g1 onto g1. g1 is changed in place
#
#   :param dict g1: genome
#   :param dict g2: genome
#   :returns: changes g1 in place"""
#   cdef list missing_chrom = set(g2.keys()) - set(g1.keys())
#   cdef int ch
#   for ch in missing_chrom:
#     g1[ch] = list(g2[ch])
#
#
# def merge_genomes(dict g1, dict g2):
#   """Given two genomes run merge_variants(c1, c2) on each chromosome. g1 has priority. In the end copy any chromosomes
#   present in g2 but not in g1 into g1 (so that g1 is a superset of g2)
#
#   :param dict g1: genome
#   :param dict g2: genome
#   :returns: new genome
#   """
#   cdef dict g3 = {}
#   for chrom, c2 in g2.iteritems():
#     g3[chrom] = merge_variants(g1.get(chrom, []), c2)
#   copy_missing_chromosomes(g3, g1)  # The previous loop merges all chr in c2.
#                                     # Now, we need to consider all chr in c1 but not in c2
#   return g3