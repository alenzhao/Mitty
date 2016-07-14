""""""
import logging

import click
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import mitty.benchmarking.lib as mlib


logger = logging.getLogger(__name__)


@click.command()
@click.version_option()
@click.argument('bdf', type=click.Path(exists=True))
@click.argument('tag', default='MQ')
@click.option('-o', default='mq.png', type=click.Path(), help='Output plot name')
@click.option('--title', help='Title of plot')
@click.option('--block-size', default=1000000, help='Block size')
@click.option('-v', count=True, help='Verbosity level')
def cli(bdf, tag, o, title, block_size, v):
  """Plot histogram of MQ (Mapping quality) or other chosen tag at given d thresholds for a given bdf.

  \b
    1. Plot frac reads incorrect vs MQ/UQ/chosen tag for different d thresholds
    2. Plot distribution of MQ/UQ/other tag values
    3. Plot mean MQ/UQ vs d_error
  """
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  H, d_r, d_be, q_be, unmapped_cnt = bin_data(bdf, tag, block_size)
  create_figure(H, d_r, d_be, q_be, unmapped_cnt, tag, title=title, show_MQ_ref=(tag=='MQ'))
  plt.savefig(o)


def bin_data(bdf_name, tag_name, block_size):
  d_max = 1000000000 + 1 # for unmapped and reads on wrong chroms
  d_range = 200
  d_bin_edges = [-1000000000.5] + np.arange(-d_range - 0.5, d_range + 1.5).tolist() + [d_max]
  qty_bin_edges = np.arange(-0.5, 71.5)
  unmapped_cnt = 0

  #                     MQ/UQ     d_range +1 is everything >= d_range
  binned_qty = np.zeros((d_range * 2 + 3, 71), dtype=int)
  st = pd.HDFStore(bdf_name, 'r')
  tn = 'bdf' if 'bdf' in st else 'crdf'
  for df in mlib.store_iter(
    bdf_name, tn, ['m1_d_error', 'm2_d_error', 'm1_mapped', 'm2_mapped', 'm1_' + tag_name, 'm2_' + tag_name], block_size):
    idx_m = (df['m1_mapped'] & 0b01 == 1)
    unmapped_cnt += idx_m.shape[0] - idx_m.sum()
    H1, _, _ = np.histogram2d(df['m1_d_error'][idx_m], df['m1_' + tag_name][idx_m], bins=(d_bin_edges, qty_bin_edges))
    idx_m = (df['m2_mapped'] & 0b01 == 1)
    unmapped_cnt += idx_m.shape[0] - idx_m.sum()
    H2, _, _ = np.histogram2d(df['m2_d_error'][idx_m], df['m2_' + tag_name][idx_m], bins=(d_bin_edges, qty_bin_edges))
    binned_qty += H1.astype(int) + H2.astype(int)

  return binned_qty, d_range, d_bin_edges, qty_bin_edges, unmapped_cnt


def create_figure(qty, d_range, d_be, qty_be, unmapped_cnt, tag_name, title='', show_MQ_ref=False):
  fig = plt.figure(1, figsize=(8, 10))
  # fig.suptitle(title)
  plt.subplots_adjust(left=0.1, hspace=0.3, top=0.95)

  ax1 = plt.subplot(411)

  mapped_read_cnt = qty.sum(axis=0)
  total_read_count = mapped_read_cnt.sum() + unmapped_cnt

  qty_vals = (qty_be[:-1] + qty_be[1:])/2
  if show_MQ_ref:
    theoretical_fw = 10.0**(-qty_vals / 10.0)
    plt.plot(qty_vals, theoretical_fw, 'k:')

  w_l = [0, 10, 100]
  fmt_l = ['k', 'b', 'g']
  label_l = ['d=0', 'd<=10', 'd<=100']

  for w, fmt, label in zip(w_l, fmt_l, label_l):
    correct_read_count = qty[d_range + 1 - w: d_range + 1 + w + 1, :].sum(axis=0)
    fw = 1.0 - correct_read_count.astype(float) / mapped_read_cnt
    # plt.plot(qty_vals, fw, fmt, drawstyle='steps-post', label=label)
    plt.plot(qty_vals, fw, fmt + '.', label=label)

  plt.legend(loc='lower left', fontsize=9)
  plt.yscale('log')
  plt.xlabel(tag_name)
  plt.ylabel('Fraction wrong')
  # plt.setp(ax1.get_xticklabels(), visible=False)
  # plt.yscale('symlog', linthreshy=0.1)

  ax2 = plt.subplot(412, sharex=ax1)
  plt.plot(qty_vals, mapped_read_cnt, 'k', drawstyle='steps-post')
  plt.text(5, total_read_count / 10,
           '{} unmapped reads ({:0.3} %)'.format(unmapped_cnt, 100.0 * unmapped_cnt/total_read_count),
           fontsize=9)
  plt.yscale('log')
  plt.xlabel(tag_name)
  plt.ylabel('Read count')

  ax3 = plt.subplot(413)
  qty_mean = np.dot(qty, qty_vals) / qty.sum(axis=1)
  plt.plot(d_be[d_range + 1 - 200: d_range + 1 + 200 + 1], qty_mean[d_range + 1 - 200: d_range + 1 + 200 + 1], 'k.')
  plt.xlabel(r'$d_{error}$')
  plt.ylabel('Mean ' + tag_name)
  plt.axvline(0, linestyle=':', color='k')
  plt.setp(ax3, xlim=[-200, 200])

  ax4 = plt.subplot(414)
  qty_tot = qty.sum(axis=1)
  plt.plot(d_be[d_range + 1 - 200: d_range + 1 + 200 + 1], qty_tot[d_range + 1 - 200: d_range + 1 + 200 + 1], 'k', drawstyle='steps-post')
  plt.text(-180, total_read_count / 10,
           '{} unmapped reads ({:0.3} %)'.format(unmapped_cnt, 100.0 * unmapped_cnt/total_read_count),
           fontsize=9)
  plt.yscale('log')
  plt.xlabel(r'$d_{error}$')
  plt.ylabel('Read count')
  plt.axvline(0, linestyle=':', color='k')
  plt.setp(ax4, xlim=[-200, 200])


if __name__ == '__main__':
  cli()