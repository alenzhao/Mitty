import numpy as np
import h5py

from mitty.version import __version__


class Population:
  """This class abstracts the storage and retrieval of the master list and samples of a population

  genome_metadata = [
    {'seq_id': 'chr1', 'seq_len': 100, 'seq_md5': 'deadbeef'},
    {'seq_id': 'chr2', 'seq_len': 200, 'seq_md5': '1337'},
  ]
  pop = Population(genome_metadata)

  The master list and samples for each chromosome are added separately. This
  """
  def __init__(self, fname=None, genome_metadata=None):
    """Load a population from file, or create a new file. Over write or store the passed master list and/or samples

    :param fname: name of the file to store/load data from. If None an in-memory file is created
    :param genome_metadata: [{seq_id, seq_len, seq_md5} ...] in same order as seen in fa.gz file
                             same format as returned by Fasta.get_seq_metadata

    The behavior is as follows:
    1. If fname is None, create a HDf5 file in memory
    2. If fname exists, open the file and load the master list
    3. If a master list is given, overwrite the existing master list if any and erase the samples
    """
    if fname is None:
      self.fp = h5py.File(name=hex(id(self)), driver='core', backing_store=False)  # Create it in memory
    else:
      self.fp = h5py.File(name=fname)
    if len(self.get_genome_metadata()) == 0:  # This is an indication that we are creating a new file
      if genome_metadata is None:
        raise RuntimeError('Creating a new Population object requires genome metadata')
      self.set_genome_metadata(genome_metadata)
      self.fp.attrs['Mitty version'] = __version__

  @staticmethod
  def get_chrom_key(chrom):
    return '/chrom_{:d}'.format(chrom)

  @staticmethod
  def get_ml_key(chrom):
    return '/chrom_{:d}/master_list'.format(chrom)

  @staticmethod
  def get_sample_key(chrom):
    return '/chrom_{:d}/samples/'.format(chrom)

  def set_genome_metadata(self, genome_metadata):
    """Save chromosome sequence metadata

    :param genome_metadata: [{seq_id, seq_len, seq_md5} ...] in same order as seen in fa.gz file
                           same format as returned by Fasta.get_seq_metadata
    """
    for n, meta in enumerate(genome_metadata):
      chrom = n + 1  # By convention we number chromosomes 1, 2, 3 ... in the same order as in the fa.gz file
      chrom_key = self.get_chrom_key(chrom)
      chrom_grp = self.fp[chrom_key] if chrom_key in self.fp else self.fp.create_group(chrom_key)
      for k, v in meta.iteritems():
        chrom_grp.attrs[k] = v

  def get_genome_metadata(self):
    """Get chromosome metadata

    :returns [{seq_id, seq_len, seq_md5} ...] same format as returned by Fasta.get_seq_metadata
    """
    return [dict(self.fp[self.get_chrom_key(n)].attrs) for n in self.get_chromosome_list()]

  def get_chromosome_metadata(self, chrom):
    return dict(self.fp[self.get_chrom_key(chrom)].attrs)

  def get_chromosome_list(self):
    return [int(ch[6:]) for ch in self.fp.keys() if ch.startswith('chrom_')]

  def set_master_list(self, chrom, master_list):
    """Replace any existing master list with this one. Erase any existing samples

    :param master_list:
    """
    assert master_list.sorted, 'Master list has not been sorted. Please check your program'
    assert len(master_list) <= 1073741823, 'Master list has more than 2^30-1 variants.'  # I want whoever gets here to mail me: kaushik.ghose@sbgenomics.com

    key = self.get_ml_key(chrom)
    if key in self.fp:
      del self.fp[key]
    sample_key = self.get_sample_key(chrom)
    if sample_key in self.fp:
      del self.fp[sample_key]
    dt = h5py.special_dtype(vlen=bytes)
    self.fp.create_dataset(self.get_ml_key(chrom), shape=master_list.variants.shape,
                           dtype=[('pos', 'i4'), ('stop', 'i4'), ('ref', dt), ('alt', dt), ('p', 'f2')],
                           data=master_list.variants, chunks=True, compression='gzip')

  def add_sample_chromosome(self, chrom, sample_name, indexes):
    """Add sample. Overwrite if name matches existing sample. No check is done to ensure list is sorted

    :param chrom:  chrom number [1, 2, 3, ...]
    :param sample_name:
    :param indexes: [(chrom, gt) ...]
    """
    key = self.get_sample_key(chrom) + '/' + sample_name
    if key in self.fp:
      del self.fp[key]
    self.fp.create_dataset(name=key, shape=indexes.shape, dtype=[('index', 'i4'), ('gt', 'i1')], data=indexes, chunks=True, compression='gzip')
    if sample_name not in self.fp.attrs.get('sample_names', []):
      self.fp.attrs['sample_names'] = list(self.fp.attrs.get('sample_names', [])) + [sample_name]

  def get_master_list(self, chrom):
    """This function loads the whole data set into memory. We have no need for chunked access right now"""
    ml = VariantList()
    if self.get_ml_key(chrom) in self.fp:
      ml.variants = self.fp[self.get_ml_key(chrom)][:]
    return ml

  def get_sample_chromosome(self, chrom, sample_name):
    """This function loads the whole data set into memory. We have no need for chunked access right now"""
    sample_key = self.get_sample_key(chrom) + '/' + sample_name
    return self.fp[sample_key][:] if sample_key in self.fp else np.array([], dtype=[('index', 'i4'), ('gt', 'i1')])

  def get_sample_names(self):
    """Return a list of sample names"""
    return self.fp.attrs['sample_names']

  def get_version(self):
    return self.fp.attrs['Mitty version']


def l2ca(l):
  """Convenience function that converts a Python list of tuples into an numpy structured array corresponding to a
  chromosome index array"""
  return np.array(l, dtype=[('index', 'i4'), ('gt', 'i1')])


class VariantList:
  """Use a numpy recarray to store a list of variants."""
  def __init__(self, pos_a=[], stop_a=[], ref_a=[], alt_a=[], p_a=[]):
    """Initialize the array from individual arrays/lists

    :param pos_a: position vector
    :param stop_a: stop vector
    :param ref_a: reference bases
    :param alt_a: alt bases
    :param p_a: probability value for the variants
    :return: Object
    """
    self.variants = np.core.records.fromarrays([pos_a, stop_a, ref_a, alt_a, p_a],
                                               dtype=[('pos', 'i4'), ('stop', 'i4'), ('ref', 'object'), ('alt', 'object'), ('p', 'f2')])
    self.sorted = False
    self.site_freq_spectrum = None

  def __len__(self):
    return self.variants.shape[0]

  def add(self, pos_a=[], stop_a=[], ref_a=[], alt_a=[], p_a=[]):
    """Add more variants to the list

    :param pos_a: position vector
    :param stop_a: stop vector
    :param ref_a: reference bases
    :param alt_a: alt bases
    :param p_a: probability value for the variants
    :return: Object
    """
    new_variants = np.core.records.fromarrays([pos_a, stop_a, ref_a, alt_a, p_a],
                                              dtype=[('pos', 'i4'), ('stop', 'i4'), ('ref', 'object'), ('alt', 'object'), ('p', 'f2')])
    self.variants = np.concatenate((self.variants, new_variants))
    self.sorted = False

  def sort(self):
    """Sort us in place by the position."""
    idx = self.variants['pos'].argsort()  # recarray sort uses all fields and is wasteful
    self.variants = self.variants[idx]
    self.sorted = True

  def balance_probabilities(self, p, f):
    """Use the ideal site probability spectrum to rescale the probability values
    :param p: probability values
    :param f: proportion

    sum(f) = 1.0 for this to work"""
    assert len(p) == len(f)
    assert abs(1.0 - sum(f)) < 1e-6
    idx = self.variants['p'].argsort()  # We need the data sorted by probability value for this to work
    n_max = len(self)
    n = 0
    for p_i, f_i in zip(p, f):
      self.variants['p'][idx[n:n + int(f_i * n_max + .5)]] = p_i  # Over index is handled gracefully
      n += int(f_i * n_max + .5)
    self.site_freq_spectrum = (p, f)

  def select(self, rng):
    """Use the rng to select variants from the master list based on their probabilities
    :param rng: a random number generator with the .rand method returning uniform random numbers (0.0, 1.0)
    :return: idx: A list of indexes into the variants indicating which have been chosen
    """
    r = rng.rand(self.variants.shape[0], 2)
    return [(r[:, 0] < self.variants['p']).nonzero()[0], (r[:, 1] < self.variants['p']).nonzero()[0]]

  def zip_up_chromosome(self, idx0, idx1):
    """Given two chromosomes, go through each copy, variant by variant, making sure they don't clash and merging any
    homozygous ones. Return us a chromosome array. Will sort master list if not sorted
    :param idx0: proposed indexes for chrom copy 0
    :param idx1: proposed indexes for chrom copy 1
    :returns: chrom, a list of tuples (index, genotype)
    """
    if not self.sorted:
      self.sort()

    # Pass 1: get rid of the colliding variants
    pos, stop = self.variants['pos'], self.variants['stop']
    z0, z1 = avoid_collisions(pos, stop, idx0), avoid_collisions(pos, stop, idx1)

    # Pass 2: merge homozygous where needed
    return merge_homozygous(pos, z0, z1)

  def generate_chromosome(self, rng):
    """Convenient wrapper around select and zip_up_chromosome
    :param rng: a random number generator with the .rand method returning uniform random numbers (0.0, 1.0)
    :returns: chrom, a list of tuples (index, genotype)
    """
    return self.zip_up_chromosome(*self.select(rng))

  def __repr__(self):
    """Fun ASCII histogram!"""
    if self.variants.shape[0] == 0:
      return '<empty>'

    if self.site_freq_spectrum is not None:
      sfs_p, sfs = self.site_freq_spectrum
      ideal_cnt = [f * self.variants.shape[0] for f in sfs]
    else:
      sfs_p = np.linspace(0, self.variants['p'].max(), num=11)  # Default is to histogram in 11 bins
      ideal_cnt = [0 for _ in range(11)]

    # Now histogram the actual data
    dp = (sfs_p[1:] - sfs_p[:-1]) / 2.0 if len(sfs_p) > 1 else [0.5]
    actual_cnt, be = np.histogram(self.variants['p'], np.concatenate(([0], sfs_p[:-1] + dp, [sfs_p[-1] + dp[-1]])))

    # We plot it as a sideways bar-graph
    size_x = min(max(actual_cnt), 80)  # columns

    #Bring the data into this grid.
    scaling_factor = float(size_x) / max(max(actual_cnt), max(ideal_cnt))
    scaled_actual = [int(v * scaling_factor + 0.5) for v in actual_cnt]
    scaled_ideal = [int(v * scaling_factor + 0.5) for v in ideal_cnt]
    rep_str = ''
    for na, sc_a, sc_i, p in zip(actual_cnt, scaled_actual, scaled_ideal, sfs_p):
      rep_str += '{:1.2f} '.format(p)
      if sc_i <= sc_a:  # The | for the ideal comes before or overlaps with the last -
        rep_str += '-' * (sc_i - 1) + ('|' if sc_i else '') + '-' * (sc_a - sc_i) + ' {:d}\n'.format(na)
      else:  # The | comes beyond the last -
        rep_str += '-' * sc_a + ' ' * (sc_i - sc_a) + '| {:d}\n'.format(na)
    return rep_str


def avoid_collisions(pos, stop, idx):
  """Remove any overlapping variants from the sequence of variants indicated by idx

  :param pos:  array of start positions of master list variants
  :param stop: array of end positions of master list variants
  :param idx:  array of indexes into the variant list
  :return: an array of non-colliding indexes
  """
  n, n_max = 0, len(idx)
  z_idx = []
  while n < n_max:
    z_idx += [idx[n]]
    n2 = n + 1
    while n2 < n_max and pos[idx[n2]] <= stop[idx[n]]:
      n2 += 1  # Collision, skip
    n = n2
  return z_idx


def merge_homozygous(pos, z0, z1):
  """Create a chromosome out of a pair of variant lists.

  :param pos:  position array from master list
  :param z0:   indexes making chrom copy 0
  :param z1:   indexes making chrom copy 1
  :return: a list of tuples (index, genotype)
  """
  n_max0, n_max1 = len(z0), len(z1)
  n0, n1 = 0, 0
  chrom = []
  while n0 < n_max0 and n1 < n_max1:
    if pos[z0[n0]] < pos[z1[n1]]:
      chrom += [(z0[n0], 0)]
      n0 += 1
      continue
    if pos[z0[n0]] > pos[z1[n1]]:
      chrom += [(z1[n1], 1)]
      n1 += 1
      continue
    # When we get here, we are equal
    if n0 < n_max0 and n1 < n_max1:  # We are equal. Are we homozygous, or just a one in a million het?
      if z0[n0] == z1[n1]:  # Yes, a hom
        chrom += [(z0[n0], 2)]
      else:  # Just two weird hets
        chrom += [(z0[n0], 0)]
        chrom += [(z1[n1], 1)]
      n0 += 1
      n1 += 1  # Lets move along

  # Now zip in the remainders
  while n0 < n_max0:
    chrom += [(z0[n0], 0)]
    n0 += 1
  while n1 < n_max1:
    chrom += [(z1[n1], 1)]
    n1 += 1

  return np.array(chrom, dtype=[('index', 'i4'), ('gt', 'i1')])