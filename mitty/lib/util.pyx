import numpy
from mitty.lib import SEED_MAX
from mitty.lib.variation import HOMOZYGOUS, HET_10, HET_01


def initialize_rngs(master_seed, n_rngs=4):
  """Return n_rngs initialized from the master_seed"""
  return [numpy.random.RandomState(seed=seed)
          for seed in numpy.random.RandomState(seed=master_seed).randint(SEED_MAX, size=n_rngs)]


def place_poisson(rng, p, end_x):
  """Given a random number generator, a probability and an end point, generate poisson distributed events. For short
  end_p this may, by chance, generate fewer locations that normal"""
  if p == 0.0:
    return numpy.array([])
  est_block_size = end_x * p * 1.2
  these_locs = rng.poisson(lam=1./p, size=est_block_size).cumsum()
  return these_locs[:numpy.searchsorted(these_locs, end_x)]


def zygosity(num_vars=0, phet=0.5, het_rng=None, copy_rng=None):
  """This function determines heterozygosity of variants.

  Parameters
  ----------
  num_vars  : int
              How many variants
  phet      : float (0.0 <= phet <= 1.0)
              Probability that a mutation is going to be heterozygous
  het_rng   : object
              Random number generator for determining heterozygosity e.g. numpy.random.RandomState(seed=1)
  copy_rng  : object
              Random number generator for determining which copy the variant is put in
              e.g. numpy.random.RandomState(seed=1)

  Returns
  -------
  het_type  : int array
              #             0      1      2      3
              gt_string = ['0/0', '0/1', '1/0', '1/1']  # The types of genotypes

  Examples
  --------
  >>> zygosity(num_vars=10, het_rng=numpy.random.RandomState(seed=1), copy_rng=numpy.random.RandomState(seed=2))
  array([2, 3, 2, 1, 2, 2, 2, 2, 1, 3], dtype=uint8)
  """
  if num_vars == 0:
    return numpy.array([])
  het_type = numpy.empty((num_vars,), dtype='u1')
  het_type.fill(HOMOZYGOUS)  # Homozygous
  idx_het, = numpy.nonzero(het_rng.rand(het_type.size) < phet)  # Heterozygous locii
  het_type[idx_het] = HET_10  # On copy 1
  het_type[idx_het[numpy.nonzero(copy_rng.rand(idx_het.size) < 0.5)[0]]] = HET_01  # On copy 2
  return het_type


cdef markov_chain_sequence_gen(
    char first_letter, float ct_mat[4][5], unsigned char *seq, unsigned long int *l, unsigned long int max_len, rng):
  """sequence_gen(float t_mat[4][5], char *alphabet[4])
  :param (char) first_letter: the first letter of the sequence
  :param (float) ct_mat: 4x5 cumulative transition probability matrix (ACGTx)
  :param (char*) seq: a max_len long string allocation. Filled out with sequence
  :param (int*) l: length of sequence
  :param (int) max_len: We pinch off a sequence that is this long
  :param rng: numpy random number generator object that has rand

  Cumulative Transition Probability matrix (x represents the break state - end of string)

     A C G T x
  A  . . . . .
  C  . . . . .
  G  . . . . .
  T  . . . . .

  Each row indicates the cumulative probability of character P (row) going to character Q (column) and thus ends with 1.
  We use cumulative probability rather than actual probability to make our code faster - we just check if the random
  number we generated is less than col0, col1, col2 and so on...
  """
  cdef:
    unsigned char last_letter = 0, n
    const char *alphabet = "ACGT"
    float r
    bint keep_running = 1

  for n in range(4):
    if first_letter == alphabet[n]:
      last_letter = n
      break

  rg = rng.rand
  l[0] = 0

  while keep_running:
    r = rg()  # Slowest part
    for n in range(5):
      if r < ct_mat[last_letter][n]:
        break
    if n == 4:
      if l[0] > 1: keep_running = 0 # We can stop if we have a length 2 sequence
    else:
      seq[l[0]] = alphabet[n]
      last_letter = n
      l[0] += 1

    if l[0] == max_len: keep_running = 0


def markov_sequences(bytes seq, ins_pts, max_len, t_mat, rng):
  """markov_sequences(seq, ins_pts, max_len, t_mat, rng)
  Return us insertions at the requested positions.
  :param (str) seq: the reference sequence. Needed for first letters of insertions
  :param (iterable) ins_pts: iterable of insertion points
  :param (4x5 list) t_mat: transition matrix
  :param rng: numpy random number generator object that has rand
  :returns a list of insertion sequences
  """
  pre_alloc_str = 'N' * max_len  # Pre-allocated string
  cdef:
    unsigned char *s = seq
    float ct_mat[4][5]
    unsigned long int l
    unsigned char *pre_string = pre_alloc_str

  # Convert transition matrix to cumulative transition matrix
  for i in range(4):
    ct_mat[i][0] = t_mat[i][0]
    for j in range(1, 5):
      ct_mat[i][j] = ct_mat[i][j-1] + t_mat[i][j]

  insertions = []
  for ip in ins_pts:
    markov_chain_sequence_gen(s[ip], ct_mat, pre_string, &l, max_len, rng)
    insertions += [pre_string[:l]]
  return insertions