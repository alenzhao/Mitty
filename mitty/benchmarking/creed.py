"""Functions to categorize reads for further analysis."""
from itertools import izip
import re

import numpy

import pyximport
pyximport.install(setup_args={"include_dirs": numpy.get_include()})
from creed_cython import *
from pysam import AlignedSegment

from mitty.lib.reads import old_style_cigar

import logging
logger = logging.getLogger(__name__)


# First ensure this rewriting works - good unit tests, with longer dels/insertions
# Then break this into smaller functions.
# Then cythonize it.
def analyze_read(read, extended=False):
  """Given a read, process the qname and read properties to determine the correct (CHROM, POS, CIGAR) and determine
  what kind of alignment errors were made on it.

  :param read: a psyam AlignedSegment object
  :param extended: whether the
  :returns read_serial, chrom, cpy, ro, pos, cigar, chrom_c, pos_c, cigar_c, unmapped

  read_serial = read_serial * 10 + 0 or 1 (for mate1 or mate2 of read) for paired reads
  read_serial = read_serial for un-paired reads
  chrom, pos -> zero indexed
  chrom_c, pos_c, cigar_c -> 1 if correct, 0 otherwise
  unmapped -> 1 if true
  """
  early_exit_value = [None] * 16

  # Don't count secondary, supplementary etc reads
  if read.flag > 255:
    return early_exit_value

  ro_m, pos_m, rl_m, cigar_m = 0, 0, 0, ''  # These are the values passed in for unpaired reads
  # We should never actually fail this, unless a tool messes badly with the qname
  try:
    #  'read_serial|chrom|copy|ro|pos|rlen|cigar|ro|pos|rlen|cigar'
    if read.is_paired:
      if read.is_read1:
        rs, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m = read.qname.split('|')
      else:
        rs, chrom, cpy, ro_m, pos_m, rl_m, cigar_m, ro, pos, rl, cigar = read.qname.split('|')
      read_serial = int(rs) * 10 + (not read.is_read1)
    else:
      rs, chrom, cpy, ro, pos, rl, cigar = read.qname.split('|')[:7]
      read_serial = int(rs)
    ro, chrom, cpy, pos, rl, ro_m, pos_m, rl_m = int(ro), int(chrom), int(cpy), int(pos), int(rl), int(ro_m), int(pos_m), int(rl_m)
  except ValueError:
    msg = 'Error processing read: {}\n'.format(read)
    # msg += '\n'.join(['{}: {}'.format(d, getattr(read, d)) for d in dir(read) if not hasattr(getattr(read, d), '__call__')])
    logger.critical(msg)
    raise RuntimeError(msg)

  chrom_c, pos_c, cigar_c, unmapped, d = 1, 1, 1, 0, MAX_D_ERROR

  if not extended:
    cigar = old_style_cigar(cigar)
    cigar_m = old_style_cigar(cigar_m)

  if read.is_unmapped:
    unmapped = 1
    chrom_c, pos_c, cigar_c = 0, 0, 0  # These are wrong by definition
  else:
    if read.reference_id != chrom - 1:
      chrom_c, pos_c = 0, 0  # chrom wrong, so pos wrong too
    else:  # Analyze the correctness by checking each breakpoint
      r = AlignedSegment()  # probably inefficient
      r.cigarstring = cigar
      cigar_ops = r.cigartuples

      # Corner case, our special cigar for indicating reads inside an insertion
      # We use S or I for this
      if cigar_ops[0][0] in [1, 4] and len(cigar_ops) == 1:  # S,I
        d = min(abs(read.pos - pos), MAX_D_ERROR)
      else:  # Go through breakpoints
        correct_pos = pos
        d_list = []
        for op, cnt in cigar_ops:
          if op in [0, 7, 8]:  # M, =, X
            d_list += [abs(read.pos - correct_pos)]
            correct_pos += cnt
          elif op == 2:  # D
            correct_pos += cnt
        d = min(d_list + [MAX_D_ERROR]) if chrom_c else MAX_D_ERROR

      if d > 0: pos_c = 0

    # This is a very strict checking of the CIGAR
    # TODO: Make this more clever
    if read.cigarstring != cigar:
      cigar_c = 0

  return read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped, d


# TODO: Remove window parameter
cigar_parser = re.compile(r'(\d+)(\D)')
def check_read(read_pos, read_cigar, correct_pos, correct_cigar, window):
  """

  :param read_pos:
  :param read_cigar:
  :param correct_pos:
  :param correct_cigar:
  :param window:
  :return: cat 3 bit value: category_code fragment indicating (pos bit, cigar bits)

  * For now the CIAGR bit is always set to 11 (fully correct) if pos is correct
  """
  cat = 0b111  # all wrong
  cigar_ops = cigar_parser.findall(correct_cigar)

  # Corner case, our special cigar for indicating reads inside an insertion
  if len(cigar_ops) == 1:
    if cigar_ops[0][1] == 'S' or cigar_ops[0][1] == 'I':
      if -window <= read_pos - correct_pos <= window:
        return 0b000

  # Not comparing CIGARs right now. Will do for the future
  for cnt, op in cigar_ops:
    if op == '=' or op == 'M' or op == 'X':
      if -window <= read_pos - correct_pos <= window:
        cat = 0b000  # all correct
        break
      correct_pos += int(cnt)
    elif op == 'D':
      correct_pos += int(cnt)
  return cat