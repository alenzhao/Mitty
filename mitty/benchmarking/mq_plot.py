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

from mitty.benchmarking.creed import MAX_D_ERROR

logger = logging.getLogger(__name__)


d_error_categories = [
  ('= 0', 0),
  ('<= 10', 10),
  ('<= 100', 100)
]


@click.command()
@click.version_option()
@click.argument('perbam', type=click.Path(exists=True))
@click.option('-o', type=click.Path(), help='Output plot name')
@click.option('--title', help='Title of plot')
@click.option('-v', count=True, help='Verbosity level')
@click.option('-p', is_flag=True, help='Show progress bar')
def cli(perbam, o, title, v, p):
  """Plot mapping quality as a function of read accuracy.

  \b
    1. Plot a binned scatter plot of read error (bp from correct pos, x-axis), mapping quality (y-axis)
    2. Plot marginals for MQ and d_error
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
  binned_mq = np.zeros((256, MAX_D_ERROR + 1), dtype=int)
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

  ax2d_pos = [left, bottom, width, height]
  rect_histx = [left, bottom_h, width, 0.2]
  rect_histy = [left_h, bottom, 0.2, height]

  # start with a rectangular Figure
  plt.figure(1, figsize=(8, 8))

  ax2d = plt.axes(ax2d_pos)
  ax_margin_MQ = plt.axes(rect_histx)
  ax_margin_d = plt.axes(rect_histy)

  ax_margin_MQ.xaxis.set_major_formatter(nullfmt)
  ax_margin_d.yaxis.set_major_formatter(nullfmt)

  mq_values = np.arange(256)
  mq_lim = (-3, binned_mq.sum(axis=1)[:-1].nonzero()[0][-1] + 3)

  d_values = np.arange(MAX_D_ERROR + 1)
  d_lim = (-10, MAX_D_ERROR + 10)

  ax2d.pcolor(mq_values, d_values, binned_mq.T,
              norm=LogNorm(vmin=1, vmax=binned_mq.max()),
              cmap='gray_r')
  plt.setp(ax2d, xlim=mq_lim, ylim=d_lim, xlabel='MQ', ylabel='|d_error| (bp)')

  # colors = [(v, v, v) for v in [0, 0.3, 0.6]]
  # for d_error_cat, color in zip(d_error_categories, colors):
  #   ax_margin_MQ.semilogy(mq_values, binned_mq[:, :d_error_cat[1] + 1].sum(axis=1), color=color, lw=2, alpha=0.71)
  ax_margin_MQ.semilogy(mq_values, binned_mq[:, :1].sum(axis=1) / float(binned_mq[:, :1].sum()), color='g', lw=3, alpha=0.71, label='|d| = 0')
  ax_margin_MQ.semilogy(mq_values, binned_mq[:, 1:].sum(axis=1) / float(binned_mq[:, 1:].sum()), color='r', lw=1, alpha=0.71, label='|d| > 0')
  plt.setp(ax_margin_MQ, xlim=mq_lim, ylim=[0, 10])
  ax_margin_MQ.legend(loc='upper center', fontsize=7)

  try:
    ax_margin_d.semilogx(binned_mq[30:, :].sum(axis=0) / float(binned_mq[30:, :].sum()), d_values, label='MQ >= 30', color='g', lw=3)
  except ValueError:
    pass
  try:
    ax_margin_d.semilogx(binned_mq[:30, :].sum(axis=0) / float(binned_mq[:30, :].sum()), d_values, label='MQ < 30', color='r', lw=1)
  except ValueError:
    pass
  plt.setp(ax_margin_d, ylim=d_lim)
  ax_margin_d.legend(loc='upper right', fontsize=7)

  # fig = plt.figure()
  # ax = fig.add_subplot(211)
  # plt.pcolor(np.arange(255), np.arange(MAX_POS + 1), binned_mq.T,
  #            norm=LogNorm(vmin=1, vmax=binned_mq.max()),
  #            cmap='gray_r')
  # plt.setp(ax, xlim=(-10, 255), ylim=(-10, MAX_POS), xlabel=' MQ', ylabel='|d_error| (bp)')
