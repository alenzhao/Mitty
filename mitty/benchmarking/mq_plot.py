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

  plot_mq(bin_mq(perbam, p), title=title)
  plt.savefig(o)


def bin_mq(perbam, p):
  """Main processing function that goes through the bam file, analyzing read alignment and writing out

  :param bam_in_fp:  Pointer to original BAM
  :param progress_bar_update_interval: how many reads to process before yielding (to update progress bar as needed)
  :return: number of reads processed
  """
  t0 = time.time()
  bam_in_fp = pysam.AlignmentFile(perbam, 'rb')
  binned_mq = np.zeros((256, 2 * MAX_D_ERROR + 2), dtype=int)
  total_read_count = bam_in_fp.mapped + bam_in_fp.unmapped  # Sadly, this is only approximate
  progress_bar_update_interval = int(0.01 * total_read_count)
  with click.progressbar(length=total_read_count, label='Reading BAM',
                         file=None if p else io.BytesIO()) as bar:
    n0 = progress_bar_update_interval
    tot_read_cnt = 0
    for tot_read_cnt, read in enumerate(bam_in_fp):
      binned_mq[read.get_tag('XM'), MAX_D_ERROR + read.get_tag('Xd')] += 1
      n0 -= 1
      if n0 == 0:
        bar.update(progress_bar_update_interval)
        n0 = progress_bar_update_interval
  t1 = time.time()
  logger.debug('Processed {:d} reads in {:2.2f}s.'.format(tot_read_cnt, t1 - t0))
  return binned_mq


def plot_mq(binned_mq, title='MQ'):
  nullfmt = NullFormatter()

  # left, width = 0.1, 0.65
  # bottom, height = 0.1, 0.65
  # bottom_h = left_h = left + width + 0.02
  #
  # ax2d_pos = [left, bottom, width, height]
  # read_cnt_hist_pos = [left + 0.2, bottom + 0.5, width - 0.4, height * 0.3]
  # rect_histx = [left, bottom_h, width, 0.2]
  # rect_histy = [left_h, bottom, 0.2, height]

  l1, b1, w1, h1 = 0.1, 0.075, 0.64, 0.48
  l2, w2 = l1 + w1 + 0.02, 0.2
  b2 = b1 + h1 + 0.02
  h2, h3 = 0.15, 0.2
  b3 = b2 + h2 + 0.02

  ax2d_pos = [l1, b1, w1, h1]
  read_cnt_hist_pos = [l1, b2, w1, h2]
  frac_correct_vs_MQ_pos = [l1, b3, w1, h3]
  d_vs_mean_MQ_pos = [l2, b1, w2, h1]

  fig = plt.figure(1, figsize=(8, 8))
  fig.suptitle(title)

  ax2d = plt.axes(ax2d_pos)
  ax_read_cnt_hist = plt.axes(read_cnt_hist_pos)
  ax_frac_correct_vs_MQ = plt.axes(frac_correct_vs_MQ_pos)
  ax_d_vs_mean_MQ = plt.axes(d_vs_mean_MQ_pos)

  ax_read_cnt_hist.xaxis.set_major_formatter(nullfmt)
  ax_frac_correct_vs_MQ.xaxis.set_major_formatter(nullfmt)
  ax_d_vs_mean_MQ.yaxis.set_major_formatter(nullfmt)

  mq_values = np.arange(256)
  mq_lim = (-3, binned_mq.sum(axis=1)[:-1].nonzero()[0][-1] + 3)

  d_values = np.arange(-MAX_D_ERROR, MAX_D_ERROR + 1)
  d_lim = (-MAX_D_ERROR - 10, MAX_D_ERROR + 10)

  ax2d.pcolor(mq_values, d_values, binned_mq[:,:-1].T,
              norm=LogNorm(vmin=1, vmax=binned_mq.max()),
              cmap='gray_r')
  plt.setp(ax2d, xlim=mq_lim, ylim=d_lim, xlabel='MQ', ylabel=r'$d_{error}$(bp)')

  ax_frac_correct_vs_MQ.step(mq_values, binned_mq[:, MAX_D_ERROR] / binned_mq.sum(axis=1).astype(float), color='k', label='d=0')
  ax_frac_correct_vs_MQ.step(mq_values, binned_mq[:, MAX_D_ERROR - 100:MAX_D_ERROR + 100].sum(axis=1) / binned_mq.sum(axis=1).astype(float), color='r', label='|d| < 100')

  ax_frac_correct_vs_MQ.plot(mq_values, 1 - 10.0**(-mq_values/10.0), color='k', linestyle=':', label=r'$10^\frac{-MQ}{10}$')
  ax_frac_correct_vs_MQ.legend(loc='lower right')
  plt.setp(ax_frac_correct_vs_MQ, xlim=mq_lim, ylim=[-0.1, 1.1], ylabel='Frac. reads\ncorrectly aligned')
  ax_frac_correct_vs_MQ.yaxis.tick_right()

  cum_rds = binned_mq[:,:-1].sum(axis=1).cumsum()
  ax_read_cnt_hist.step(mq_values, cum_rds, color='gray', lw=2)
  rds = binned_mq[:,:-1].sum(axis=1)
  ax_read_cnt_hist.step(mq_values, rds, color='k')
  plt.setp(ax_read_cnt_hist, xlim=mq_lim, ylim=(rds.min()/5., 5 * cum_rds.max()), ylabel='Read counts', yscale='log')
  ax_read_cnt_hist.yaxis.tick_right()
  ax_read_cnt_hist.yaxis.grid(True)

  ax_d_vs_mean_MQ.plot(np.dot(mq_values, binned_mq[:, :-1])/binned_mq[:, :-1].sum(axis=0).astype(float), d_values, 'k.')
  plt.setp(ax_d_vs_mean_MQ, ylim=d_lim, xlim=mq_lim, xlabel='Mean MQ')
