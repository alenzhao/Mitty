import string

import numpy

from mitty.lib import SEED_MAX


DNA_complement = string.maketrans('ATCGN', 'TAGCN')


def initialize_rngs(unsigned long master_seed, int n_rngs=4):
  """Return n_rngs initialized from the master_seed"""
  return [numpy.random.RandomState(seed=seed)
          for seed in numpy.random.RandomState(seed=master_seed).randint(SEED_MAX, size=n_rngs)]


def place_poisson_seq(rng, float p, unsigned long start_x, unsigned long end_x, bytes seq):
  """Given a random number generator, a probability and an end point, generate poisson distributed events. Skip bases
  that don't belong  For short end_p this may, by chance, generate fewer locations that normal"""
  if p == 0.0:
    return numpy.array([])

  cdef:
    char *s = seq
    unsigned long est_block_size = <unsigned long>(<float>end_x * p * 1.2)
    unsigned long idx

  these_locs = rng.geometric(p=p, size=est_block_size).cumsum()
  return numpy.array([idx for idx in these_locs[numpy.searchsorted(these_locs, start_x):numpy.searchsorted(these_locs, end_x)] if s[idx] != 'N'], dtype='i4')


cdef unsigned char sub_base(unsigned char orig_base, unsigned char sub_mat[85][3], float ct_mat[85][3], float r):
  # sub_matrix
  #          0  1  2
  # A (65)   C  G  T
  # C (67)   A  G  T
  # G (71)   A  C  T
  # T (84)   A  C  G

  cdef:
    unsigned char i
  for i in range(3):
    if r < ct_mat[orig_base][i]:
      break
  return sub_mat[orig_base][i]


def base_subs(bytes seq, sub_pts, t_mat, rng):
  """
  t_mat ->
     A C G T
  A  . x x x
  C  x . x x
  G  x x . x
  T  x x x .
  """
  cdef:
    char *s = seq
    unsigned char i, j, ob, sb
    unsigned char sub_mat[85][3]
    float ct_mat[85][3]

  # Create the ct_mat from the basic t_mat
  for i, ob in enumerate([65, 67, 71, 84]):
    jj = 0
    for j, sb in enumerate([65, 67, 71, 84]):
      if i == j: continue
      ct_mat[ob][jj] = t_mat[i][j]  # First pass, compact the matrix
      sub_mat[ob][jj] = sb  #
      jj += 1
    # Second pass, do cumulative probabilities
    ct_mat[ob][1] += ct_mat[ob][0]
    ct_mat[ob][2] += ct_mat[ob][1]

  r = rng.rand(len(sub_pts))
  cdef unsigned long q
  return [<bytes> sub_base(s[sub_pts[q]], sub_mat, ct_mat, r[q]) for q in range(len(sub_pts))]


def add_p_end_to_t_mat(t_mat, p_end):
  """Given p_end, incorporate it into the t_mat."""
  p_end_1 = 1.0 - p_end
  return [[p_end_1 * (t_mat[i][j]/sum(t_mat[i])) for j in range(4)] + [p_end] for i in range(4)]


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
  :param (4x5 list) t_mat: transition matrix, including prob of termination
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