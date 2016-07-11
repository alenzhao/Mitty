import json
import sys
import os
import time
import base64
import logging
from multiprocessing import Process, Queue, current_process, freeze_support
import subprocess

import pysam
import click

from mitty.version import __version__


logger = logging.getLogger(__name__)

__process_stop_code__ = 'SETECASTRONOMY'


@click.command()
@click.version_option()
@click.argument('fasta_ann', type=click.Path(exists=True))
@click.argument('fastq', type=click.Path(exists=True))
@click.argument('bam', type=click.Path())
@click.option('-t', default=2, help='Number of threads')
@click.option('--fastq2', type=click.Path(exists=True), help='File2 of a pair of FASTQs')
@click.option('--interleaved-fq', is_flag=True, default=True, help='Is the FASTQ an interleaved. Ignored if fastq2 is supplied')
@click.option('--sample-name', help="Name of sample", default='S')
@click.option('--block-size', default=1000000, help='Maximum number of templates held in memory')
@click.option('--max-templates', type=int, help='Maximum number of templates to process (For debugging)')
@click.option('-v', count=True, help='Verbosity level')
def cli(fasta_ann, fastq, bam, t, fastq2, sample_name, interleaved_fq, block_size, max_templates, v):
  """Given a FASTA.ann file and FASTQ made of simulated reads,
  construct a perfectly aligned "god" BAM from them.
  This is useful for testing variant callers.

  The program uses the fasta.ann file to construct the BAM header
  """
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  if not interleaved_fq:
    raise RuntimeError('Sorry, non-interleaved FASTQ reading not implemented yet')

  rg_id = base64.b64encode(' '.join(sys.argv))
  bam_hdr = construct_header(fasta_ann, rg_id=rg_id, sample=sample_name)

  process_multi_threaded(fastq, bam, bam_hdr, t, max_templates)
  # process_single_threaded(fastq, bam, bam_hdr, max_templates)  (Only for debugging)


def construct_header(fasta_ann, rg_id, sample='S'):
  return {
    'HD': {'VN': '1.0'},
    'PG': [{'CL': ' '.join(sys.argv),
            'ID': 'mitty-god-aligner',
            'PN': 'god-aligner',
            'VN': __version__}],
    'RG': [{'ID': rg_id, 'SM': sample}],
    'SQ': parse_ann(fasta_ann)
  }


def parse_ann(fn):
  """Given a fasta.ann file name parse it."""
  logger.debug('Parsing {} for sequence header information'.format(fn))
  ln = open(fn, 'r').readlines()[1:]  # Don't need the first line
  return [
    {
      'SN': ln[n].split()[1],
      'LN': int(ln[n + 1].split()[1])
    }
  for n in range(0, len(ln), 2)]


def process_multi_threaded(fastq, bam, bam_hdr, t, max_templates):
  in_queue = Queue()

  # Start worker processes
  logger.debug('Starting {} threads'.format(t))
  file_fragments = ['frag-{:03}-{}'.format(i, bam) for i in range(t)]
  p_list = [Process(target=process_worker, args=(file_fragments[i], bam_hdr, in_queue)) for i in range(t)]
  for p in p_list:
    p.start()

  # Burn through file
  logger.debug('Starting to read FASTQ file')
  read_fastq(fastq, True, in_queue, max_templates=max_templates)

  # Tell child processes to stop
  logger.debug('Stopping child processes')
  for i in range(t):
    in_queue.put(__process_stop_code__)

  # Wait for them to finish
  for p in p_list:
    p.join()

  # Get samtools to combine the bams together
  samtools_cat(bam, file_fragments)


def process_single_threaded(fastq, bam, bam_hdr, max_templates):
  in_queue = Queue()

  # Burn through file
  logger.debug('Starting to read FASTQ file')
  read_fastq(fastq, True, in_queue, max_templates=max_templates)
  # Place the magic stop code
  in_queue.put(__process_stop_code__)

  # Process reads
  process_worker(bam, bam_hdr, in_queue)


def read_fastq(fastq, paired_end, in_queue, max_templates=None):
  """

  :param fastq:
  :param input:  object where data can put pushed to using put
  :param max_templates:
  :return:
  """
  f_size = os.stat(fastq).st_size
  template_notification_interval = 100000
  template_counter = 0
  t0 = time.time()
  with open(fastq, 'r') as fp:
    template = []
    read = []
    reads_in_template_count = 2 if paired_end else 1
    line_cnt = 4
    for ln in fp:
      read.append(ln[:-1])
      line_cnt -= 1
      if line_cnt == 0:
        template.append(read)
        read = []
        line_cnt = 4
        reads_in_template_count -= 1
        if reads_in_template_count == 0:
          in_queue.put(template)
          template_counter += 1
          if template_counter % template_notification_interval == 0:
            logger.debug('Processed {} templates ({:0.7} templates/s, {:0.5}/{:0.5} MB)'.format(template_counter, template_counter/(time.time() - t0), fp.tell()/1e6, f_size/1e6))
          if max_templates is not None and template_counter >= max_templates:
            logger.debug('Stopping at {} templates, as asked for'.format(max_templates))
            break
          template = []
          reads_in_template_count = 2 if paired_end else 1

  t1 = time.time()
  logger.debug('Took {:0.5}s to process {} templates ({:0.7} tmpl/s)'.format(t1 - t0, template_counter, template_counter/(t1 - t0)))


def process_worker(bam_fname, bam_hdr, in_queue):
  """Create a bam fragment. This is designed to be a worker process for a multiprocessing
  pool, but can be tested without recourse to multiprocessing

  :param bam_fname:  output file name the reads will be written to
  :param bam_hdr:    the bam header (as a dict)
  :param in_queue:      an object with a get method that returns templates when called with next()
  :return:
  """
  logger.debug('Writing to {}'.format(bam_fname))
  fp = pysam.AlignmentFile(bam_fname, 'wb', header=bam_hdr)
  for template in iter(in_queue.get, __process_stop_code__):
    process_template(template, fp)
  fp.close()
  logger.debug('Shutting down thread for {}'.format(bam_fname))


def process_template(template, fp):
  """Given a template, write out the perfect alignments to file

  :param template: [x1, x2, ... ] where xi is a group of four lines from the FASTQ file
                   and, e.g. [x1, x2] constitute a pair, if the input is paired end
  :param fp: pysam file pointer to write out
  :return:
  """
  paired_end = len(template) == 2
  qname = template[0][0][1:]  # First line is qname, startswith an @
  if paired_end:
    qname = qname[:-2]  # Paired-end qnames end with /1, /2 etc.
    rs, chrom_s, cpy_s, ro_s1, pos_s1, rl_s1, cigar1, ro_s2, pos_s2, rl_s2, cigar2 = qname.split('|')
  else:
    rs, chrom_s, cpy_s, ro_s1, pos_s1, rl_s1, cigar1 = qname.split('|')

  r1 = pysam.AlignedSegment()
  r1.qname = qname

  r1.reference_id = int(chrom_s) - 1
  r1.pos = int(pos_s1)
  r1.cigarstring = cigar1
  r1.seq = template[0][1]  # Second line is the base sequence
  r1.qual = template[0][3]  # Fourth line is the base quality
  r1.is_paired = False

  # TODO: Fix this in the simulator
  cigar_ops = r1.cigartuples
  # Corner case, our special cigar for indicating reads inside an insertion
  # We use S or I for this
  if cigar_ops[0][0] in [1, 4] and len(cigar_ops) == 1:  # S,I
    r1.cigarstring = '{}I'.format(r1.rlen)
    r1.is_unmapped = 1

  if paired_end:
    r1.is_paired = True
    r1.is_read1 = True

    r2 = pysam.AlignedSegment()
    r2.qname = qname
    r2.reference_id = int(chrom_s) - 1
    r2.pos = int(pos_s2)
    r2.cigarstring = cigar2
    r2.seq = template[1][1]  # Second line is the base sequence
    r2.qual = template[1][3]  # Fourth line is the base quality

    # TODO: Fix this in the simulator
    cigar_ops = r2.cigartuples
    # Corner case, our special cigar for indicating reads inside an insertion
    # We use S or I for this
    if cigar_ops[0][0] in [1, 4] and len(cigar_ops) == 1:  # S,I
      r2.cigarstring = '{}I'.format(r2.rlen)

    r2.is_paired = True
    r2.is_read1 = False
    r2.is_read2 = True

    r1.pnext = r2.pos
    r2.pnext = r1.pos

    r1.rnext = r2.reference_id
    r2.rnext = r1.reference_id

    fp.write(r1)
    fp.write(r2)
  else:
    fp.write(r1)


def samtools_cat(bam, file_fragments):
  cmd = ['samtools', 'cat', '-o', bam] + file_fragments
  subprocess.call(cmd)
  # Now cleanup the fragments
  # for frag in file_fragments:
  #   os.remove(frag)


if __name__ == "__main__":
  cli()