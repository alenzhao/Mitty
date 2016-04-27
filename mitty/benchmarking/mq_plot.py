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

  ax_margin_MQ.plot(mq_values, binned_mq[:, :1].sum(axis=1) / binned_mq.sum(axis=1).astype(float), color='k', label='data')
  ax_margin_MQ.plot(mq_values, 1 - 10.0**(-mq_values/10.0), color='k', linestyle=':', label=r'$10^\frac{-MQ}{10}$')
  ax_margin_MQ.legend(loc='lower center')
  plt.setp(ax_margin_MQ, xlim=mq_lim, ylim=[-0.1, 1.1], title='Frac. reads correctly aligned')

  ax_margin_d.plot(np.dot(mq_values, binned_mq)/binned_mq.sum(axis=0).astype(float), d_values, 'ko')
  plt.setp(ax_margin_d, ylim=d_lim, xlim=mq_lim, title='Mean MQ')
