"""Basic plugin producing single ended reads

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
__example_param_text = """
{
  'read_len': 100,          # length of each read
  'max_p_error': 0.01,      # Maximum error rate at tip of read
  'k': 20,                  # exponent
}
"""

_description = __doc__ + '\nExample parameter set:\n' + __example_param_text

_example_params = eval(__example_param_text)

from itertools import izip

import numpy as np

import mitty.lib.util as mutil
from mitty.plugins.reads.base_plugin import ReadModel

import logging
logger = logging.getLogger(__name__)


class Model(ReadModel):
  """Stock read plugin that approximates Illumina reads"""
  def __init__(self, read_len=100, max_p_error=0.01, k=20):
    """Initialization

    :param read_len:
    :param max_p_error:
    :param k:
    :return:
    """
    self.read_len = read_len
    l = np.linspace(1e-10, 1, self.read_len)
    self.error_profile = max_p_error * (np.exp(k*l) - 1)/(np.exp(k) - 1)
    self.phred = ''.join([chr(int(33 + max(0, min(-10*np.log10(p), 93)))) for p in self.error_profile])
    ReadModel.__init__(self, False)

  def get_reads(self, seq, seq_c, start_base=0, end_base=None, coverage=0.01, corrupt=False, seed=1):
    """The main simulation calls this function.

    :param seq:      forward sequence
    :param seq_c:    complement of sequence
    :param start_base: base to start taking reads from
    :param end_base: base to stop
    :param coverage: coverage
    :param corrupt:  T/F whether we should compute corrupted read or not
    :param seed:     random number generator seed
    :return: reads, paired

    reads is a numpy recarray with the following fields
      'start_a'  -> start index on the seq
      'read_len' -> length of the read
      'rev_complement' ->  0 or 1, indicating if the read is from the forward or reverse strand
      'perfect_read'     -> perfect read sequence
      'corrupted_read'   -> corrupted read sequence

    paired indicates if the reads are in pairs or not
    """
    end_base = end_base or len(seq)
    p_template = coverage / float(self.read_len)
    template_cnt = int(p_template * (end_base - start_base))

    template_loc_rng, read_order_rng, error_loc_rng, base_choice_rng = mutil.initialize_rngs(seed, 4)
    template_locs = template_loc_rng.randint(start_base, end_base + 1, size=template_cnt)
    idx = ((template_locs + self.read_len < end_base) & (start_base < template_locs - self.read_len)).nonzero()[0]
    template_locs = template_locs[idx]
    read_order = read_order_rng.randint(2, size=template_locs.shape[0])  # Which read comes first?

    reads = np.recarray(dtype=ReadModel.dtype, shape=template_locs.shape[0])
    r_start, r_o, pr = reads['start_a'], reads['read_order'], reads['perfect_reads']
    reads['read_len'] = self.read_len
    reads['read_order'] = read_order
    r_len = self.read_len

    idx_fwd = (read_order == 0).nonzero()[0]
    idx_rev = read_order.nonzero()[0]

    r_start[idx_fwd] = template_locs[idx_fwd]
    r_start[idx_rev] = template_locs[idx_rev] - r_len

    for n in xrange(reads.shape[0]):
      if r_o[n] == 0:  # Forward strand
        pr[n] = seq[r_start[n]:r_start[n] + r_len]
      else:  # Reverse strand
        pr[n] = seq_c[r_start[n]:r_start[n] + r_len][::-1]

    if corrupt:
      self.corrupt_reads(reads, error_loc_rng, base_choice_rng)
    return reads, self.paired

  def corrupt_reads(self, reads, error_loc_rng, base_choice_rng):
    """Corrupt reads

    :param reads:   with the start_a, read_len and read_order fields filled out
    :return: Fill in corrupted reads in place
    """
    pr, cr, ph = reads['perfect_reads'], reads['corrupt_reads'], reads['phred']
    p = error_loc_rng.rand(reads.shape[0], self.read_len)
    idx = [np.where(p[n, :] < self.error_profile)[0] for n in xrange(reads.shape[0])]
    idx_tot = 0
    for n in xrange(reads.shape[0]):
      idx_tot += idx[n].size
    corrupt_bases = base_choice_rng.choice(['A','C','G','T'], size=idx_tot, replace=True, p=[.3, .2, .2, .3]).tostring()

    offset = 0
    for n in xrange(reads.shape[0]):
      ph[n] = self.phred
      #idx, = np.where(p[n, :] < self.error_profile)
      #corrupt_bases = base_choice_rng.choice(['A','C','G','T'], size=idx.size, replace=True, p=[.3, .2, .2, .3]).tostring()
      b = bytearray(pr[n])
      for m, i in enumerate(idx[n]):
        b[i] = corrupt_bases[offset +m]
      cr[n] = str(b)
      offset += idx[n].size


def self_test():
  """Basic self test"""
  import string
  DNA_complement = string.maketrans('ATCGN', 'TAGCN')
  seq = 'ATGTCGCCGGGCGCCATGCGTGCCGTTGTTCCCATTATCCCATTCCTTTTGGTTCTTGTCGGTGTATCGGGGGTTCCCACCAACGTCTCCTCCACCACCCAACCCCAACTCCAGACCACCGGTCGTCCCTCGCATGAAGCCCCCAACATGACCCAGACCGGCACCACCGACTCTCCCACCGCCATCAGCCTTACCACGCCCGACCACACACCCCCCATGCCAAGTATCGGACTGGAGGAGGAGGAAGAGGAGGAGGGGGCCGGGGATGGCGAACATCTTGAGGGGGGAGATGGGACCCGTGACACCCTACCCCAGTCCCCGGGTCCAGCCGTCCCGTTGGCCGGGGATGACGAGAAGGACAAACCCAACCGTCCCGTAGTCCCACCCCCCGGTCCCAACAACTCCCCCGCGCGCCCCGAGACCAGTCGACCGAAGACACCCCCCACCAGTATCGGGCCGCTGGCAACTCGACCCACGACCCAACTCCCCTCAAAGGGGCGACCCTTGGTTCCGACGCCTCAACATACCCCGCTGTTCTCGTTCCTCACTGCCTCCCCCGCCCTGGACACCCTCTTCGTCGTCAGCACCGTCATCCACACCTTATCGTTTTTGTGTATTGTTGCGATGGCGACACACCTGTGTGGCGGTTGGTCCAGACGCGGGCGACGCACACACCCTAGCGTGCGTTACGTGTGCCTGCCGCCCGAACGCGGGTAG'
  seq_c = string.translate(seq, DNA_complement)
  mdl = Model(4, 8, 2, max_p_error=1)
  rd, paired = mdl.get_reads(seq, seq_c, start_base=0, end_base=len(seq), coverage=.00001, corrupt=True)
  assert type(rd) == np.core.records.recarray  # Basically, the previous code should just run
  assert paired == False

if __name__ == "__main__":
  print _description