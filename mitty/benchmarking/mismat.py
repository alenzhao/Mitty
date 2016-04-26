"""Given a perfectbam/indelbam save a .csv file of read assignments by chromosome, including mismatch."""
import time
import io
import cPickle
import json

import click
import pysam
import numpy as np

import mitty.version
import mitty.lib.variants as vr
import mitty.benchmarking.creed as creed

import logging
logger = logging.getLogger(__name__)


@click.command()
@click.argument('inbam', type=click.Path(exists=True))
@click.argument('outcsv', type=click.Path())
@click.option('-p', is_flag=True, help='Show progress bar')
@click.option('-v', count=True, help='Verbosity level')
def cli(inbam, outcsv, p, v):
  """Given a perfectbam/indelbam save a .csv file of read assignments by chromosome, including mismatch."""
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  bam_fp = pysam.AlignmentFile(inbam, 'rb')
  n_chroms = len(bam_fp.header['SQ'])
  cnt_mat = np.zeros((n_chroms, (n_chroms + 1)), dtype=int)

  total_read_count = bam_fp.mapped + bam_fp.unmapped  # Sadly, this is only approximate
  progress_bar_update_interval = int(0.01 * total_read_count)

  cnt = 0
  t0 = time.time()
  with click.progressbar(
    length=total_read_count, label='Analyzing mapping pattern', file=None if p else io.BytesIO()) as bar:

    # TAG TYPE VALUE
    # XR  i    Aligned chromosome
    # XP  i    Aligned pos
    # Xf  i    mapped correctly/incorrectly, unmapped

    for cnt, r in enumerate(bam_fp):
      c_chrom, a_chrom = r.reference_id, r.get_tag('XR')
      cnt_mat[c_chrom, a_chrom + 1] += 1  # First column is unmapped

      if (cnt + 1) % progress_bar_update_interval == 0:
        bar.update(progress_bar_update_interval)
    bar.update(progress_bar_update_interval)
  t1 = time.time()

  logger.debug('Analyzed {:d} reads in {:2.2f}s'.format(cnt, t1 - t0))

  # Now save the output CSV. Only save those rows (chromosomes) that had at least a single originating read
  # header = [i for i in range(n_chroms)]
  # rows = [cnt_mat[r, :] for r in range(n_chroms) if cnt_mat[r, :].sum()]

  with open(outcsv, 'w') as fp:
    # Header
    fp.write(','.join(['Source Chrom'] + [str(i) if i > 0 else 'Un-mapped' for i in range(n_chroms + 1)]) + '\n')
    for r in range(n_chroms):
      if cnt_mat[r, :].sum():
        fp.write(','.join([str(r + 1)] + [str(i) for i in cnt_mat[r, :]]) + '\n')


if __name__ == '__main__':
  cli()