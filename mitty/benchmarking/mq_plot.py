import logging
import time
import io

import pysam
import numpy as np
import click
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import NullFormatter

from mitty.benchmarking.perfectbam import MAX_POS

logger = logging.getLogger(__name__)


@click.command()
@click.version_option()
@click.argument('perbam', type=click.Path(exists=True))
@click.option('-o', type=click.Path(), help='Output plot name')
@click.option('--title', help='Title of plot')
@click.option('-v', count=True, help='Verbosity level')
@click.option('-p', is_flag=True, help='Show progress bar')
def cli(perbam, o, title, v, p):
  """Plot mapping quality as a function of read accuracy.

  1. Plot a binned scatter plot of read error (bp from correct pos, x-axis), mapping quality (y-axis)
  2. Plot a marginal with MQ histogram
  """
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  plot_mq(bin_mq(perbam, p))
  plt.savefig(o)


def bin_mq(perbam, p):
  """Main processing function that goes through the bam file, analyzing read alignment and writing out

  :param bam_in_fp:  Pointer to original BAM
  :param progress_bar_update_interval: how many reads to process before yielding (to update progress bar as needed)
  :return: number of reads processed
  """
  t0 = time.time()
  bam_in_fp = pysam.AlignmentFile(perbam, 'rb')
  binned_mq = np.zeros((255, MAX_POS + 1), dtype=int)
  total_read_count = bam_in_fp.mapped + bam_in_fp.unmapped  # Sadly, this is only approximate
  progress_bar_update_interval = int(0.01 * total_read_count)
  with click.progressbar(length=total_read_count, label='Reading BAM',
                         file=None if p else io.BytesIO()) as bar:
    n0 = progress_bar_update_interval
    tot_read_cnt = 0
    for tot_read_cnt, read in enumerate(bam_in_fp):
      binned_mq[read.mapping_quality, read.get_tag('Xd')] += 1
      n0 -= 1
      if n0 == 0:
        bar.update(progress_bar_update_interval)
        n0 = progress_bar_update_interval
  t1 = time.time()
  logger.debug('Processed {:d} reads in {:2.2f}s.'.format(tot_read_cnt, t1 - t0))
  return binned_mq


def plot_mq(binned_mq):
  nullfmt = NullFormatter()

  left, width = 0.1, 0.65
  bottom, height = 0.1, 0.65
  bottom_h = left_h = left + width + 0.02

  rect_scatter = [left, bottom, width, height]
  rect_histx = [left, bottom_h, width, 0.2]
  rect_histy = [left_h, bottom, 0.2, height]

  # start with a rectangular Figure
  plt.figure(1, figsize=(8, 8))

  ax2d = plt.axes(rect_scatter)
  axHistx = plt.axes(rect_histx)
  axHisty = plt.axes(rect_histy)

  axHistx.xaxis.set_major_formatter(nullfmt)
  axHisty.yaxis.set_major_formatter(nullfmt)

  mq_values = np.arange(255)
  mq_lim = (-10, 255)

  d_values = np.arange(MAX_POS + 1)
  d_lim = (-10, MAX_POS + 10)

  ax2d.pcolor(mq_values, d_values, binned_mq.T,
              norm=LogNorm(vmin=1, vmax=binned_mq.max()),
              cmap='gray_r')
  plt.setp(ax2d, xlim=mq_lim, ylim=d_lim, xlabel='MQ', ylabel='|d_error| (bp)')

  axHistx.semilogy(mq_values, binned_mq.sum(axis=1))
  plt.setp(axHistx, xlim=mq_lim)

  axHisty.semilogx(binned_mq.sum(axis=0), d_values)
  plt.setp(axHisty, ylim=d_lim)

  # fig = plt.figure()
  # ax = fig.add_subplot(211)
  # plt.pcolor(np.arange(255), np.arange(MAX_POS + 1), binned_mq.T,
  #            norm=LogNorm(vmin=1, vmax=binned_mq.max()),
  #            cmap='gray_r')
  # plt.setp(ax, xlim=(-10, 255), ylim=(-10, MAX_POS), xlabel=' MQ', ylabel='|d_error| (bp)')
