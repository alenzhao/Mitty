"""This script takes in a fastq file generated by vcf2reads and aligns it perfectly based on the cheat answers stored
in the qnames. This is useful if we want to test a tool assuming a perfect alignment of the reads.

Commandline::

  Usage:
    reads2bam  [-p] --fa_dir=FADIR  --fastq=FASTQ  --bam=BAM  [-v|-V]

  Options:
    -p                      If set indicates fastq is interleaved paired
    --fa_dir=FADIR          Directory where genome is located
    --fasta=FASTQ           Input fasta file
    --bam=BAM               Output bam file
    -v                      Dump detailed logger messages
    -V                      Dump very detailed logger messages
"""
__version__ = '1.0.0'
import docopt
import os
import pysam
from lib.genome import FastaGenome
import logging
logger = logging.getLogger(__name__)

import string
DNA_complement = string.maketrans('ATCGN', 'TAGCN')


def sort_and_index_bam(bamfile):
  """Do the filename gymnastics required to end up with a sorted, indexed, bam file."""
  # samtools sort adds a '.bam' to the end of the file name.
  os.rename(bamfile, 'temp.bam')
  pysam.sort('temp.bam', os.path.splitext(bamfile)[0])
  pysam.index(bamfile)
  os.remove('temp.bam')


def align(in_fastq, out_bam, seq_dir = '', paired=False):
  def interpret_read_qname(qname, template_order):
    """chrom:copy|rN|D|POS1|CIGAR1|POS2|CIGAR2"""
    qn = qname.split('|')
    return int(qn[0][:-2]), qn[2], int(qn[3 if not template_order else 5]), qn[4 if not template_order else 6]

  def process_read(read, template_order):
    correct_chrom_no, dirn, correct_pos, correct_cigar = interpret_read_qname(read.name, template_order)
    a_read = pysam.AlignedRead()
    a_read.pos = correct_pos - 1  # BAM files are zero indexed
    a_read.tid = correct_chrom_no - 1
    a_read.qname = read.name
    a_read.cigarstring = correct_cigar
    if (template_order and dirn == '>') or (not template_order and dirn == '<'):
      a_read.seq = read.sequence[::-1].translate(DNA_complement)
      a_read.flag = 0x10  # 0x10 flag has to be set to indicate reverse complement
    else:
      a_read.seq = read.sequence
      a_read.flag = 0x0
    a_read.qual = read.quality
    a_read.mapq = 100  # It's better to set this
    return a_read

  def get_read(fq, _paired):
    while 1:
      if _paired:
        r1, r2 = process_read(fq.next(), 0), process_read(fq.next(), 1)
        # 0b01 flag has to be set to indicate multiple segments
        # 0b10 flag has to be set to indicate each segment mapped
        r1.flag |= 0x43  # 1st segment
        r2.flag |= 0x83  # last segment
        # If we don't set template len and pnexts Tablet doesn't show us the mate pairs properly
        r1.tlen = r2.tlen = r2.pos + r2.rlen - r1.pos
        r1.pnext, r2.pnext = r2.pos, r1.pos
        r1.rnext = r2.rnext = r1.tid
        yield [r1, r2]
      else:
        yield [process_read(fq.next(), 0)]

  ref = FastaGenome(seq_dir=seq_dir)
  bam_hdr = {'HD': {'VN': '1.4'},
             'SQ': [{'LN': seq_len, 'SN': seq_id.split(' ')[0]} for seq_id, seq_len in ref.genome_header()]}
  #  Tablet and other viewers get confused by spaces in the seq_id. Losers
  fastq = pysam.Fastqfile(in_fastq)
  out_bamfile = pysam.Samfile(out_bam, 'wb', header=bam_hdr)

  for cnt, (reads) in enumerate(get_read(fastq, paired)):
    for read in reads:
      out_bamfile.write(read)

    if not cnt % 100000:
      logger.debug('{:d} templates done'.format(cnt))

  logger.debug('Aligned {:d} templates'.format(cnt))
  out_bamfile.close()
  sort_and_index_bam(out_bam)


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)

  if args['-V']:
    logging.basicConfig(level=logging.DEBUG)
  else:
    logging.basicConfig(level=logging.WARNING)

  if args['-v']:
    logger.setLevel(logging.DEBUG)

  align(args['--fastq'], args['--bam'], args['--fa_dir'], paired=True if args['-p'] else False)