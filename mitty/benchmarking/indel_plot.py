"""Plot layout:

Repeat for novel and known, separate figures
Ref info is duplicated
There are three traces on the indel plot (top row)
one for d_err = 0, d_err < 10, d_err < 100

----------------------------
|         |   |   |        |
|   DEL   | R | S |   INS  |
|         |   |   |        |
----------------------------
|   read counts            |
----------------------------
|   var counts             |
----------------------------

"""


import itertools
import cPickle

import click
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import mitty.benchmarking.creed as creed

import logging
logger = logging.getLogger(__name__)


d_error_categories = [
  ('= 0', 0),
  ('<= 10', 10),
  ('<= 100', 100)
]


@click.command()
@click.version_option()
@click.argument('indelpkl', type=click.Path(exists=True))
@click.argument('outsuffix', type=click.Path())
@click.option('--indel-range', default=100, help='Range of indels to plot')
@click.option('--win', default=0, help='Size of median filter window to smooth plots')
@click.option('--title', help='Title', default='Aligner accuracy')
def cli(indelpkl, outsuffix, indel_range, win, title):
  """pickled data -> plot:

  \b
  Given a pkl file with alignment accuracy parametrized by the sample indel lengths,
  generate a plot."""
  win += (win + 1)% 2  # Ensure win is odd
  read_counts = cPickle.load(open(indelpkl, 'r'))

  for sl in ['known', 'novel']:
    plot_fig(creed.slice_read_counts(read_counts, sl), indel_range, win, '{}-{}'.format(title, sl))
    out_name = '{}-{}'.format(sl, outsuffix)
    plt.savefig(out_name)
    print('Saved {} variant plot in {}'.format(sl, out_name))


def plot_fig(read_counts, indel_range, win, title):
  fig, ax = setup_figure()
  axes_spec = build_axes_specs(read_counts, indel_range, title)
  colors = [(v, v, v) for v in [0, 0.3, 0.6, 0.8]]
  for d_error_cat, color in zip(d_error_categories, colors):
    plot_alignment_accuracy(read_counts, d_error_cat, color=color, ax=ax)
  plot_read_count(read_counts, color='k', lw=2, ax=ax)
  plot_var_count(read_counts, color='k', lw=2, ax=ax)
  decorate_axes(ax, axes_spec)


def setup_figure():
  fig = plt.figure(figsize=(2 * 3.385, 2 * 3))  # two column figure for bio-informatics
  plt.subplots_adjust(left=0.15, bottom=0.1, right=0.98, top=0.93, wspace=0.05, hspace=0.01)
  gs = plt.GridSpec(3, 3,
                    width_ratios=[6, 0.9, 6],
                    height_ratios=[3.5, 1, 1])
  ax = {k: plt.subplot(g) for k, g in
        zip([''.join(e) for e in itertools.product(['A', 'B', 'C'], ['-DEL', '-REF/SNP', '-INS'])], gs)}
  return fig, ax


def build_axes_specs(read_counts, indel_range, title):
  rc = read_counts

  rd_cnt = np.concatenate([rc['v_r_cnt'].sum(axis=1), [rc['ref_r_cnt'].sum()]])
  max_read_cnt, min_read_cnt = max(2, rd_cnt.max()), min(1, rd_cnt.min())
  max_indel_cnt, min_indel_cnt = max(2, rc['v_cnt'].max()), min(1, rc['v_cnt'].min())
  indel_range = indel_range or rc['v_len'].max()
  axes_specs = {
    'title': title,
    'x-lim': [-indel_range, indel_range],
    'x-ticks': [-indel_range, -1, 0, 1, indel_range],
    'A': {'y-lim': [0, 105], 'y-ticks': [25, 50, 75, 95, 100]},
    'B': {'y-lim': [min_read_cnt / 2, max_read_cnt * 2],
          'y-ticks': [min_read_cnt, max_read_cnt],
          'y-t-labels': [min_read_cnt, max_read_cnt]},
    'C': {'y-lim': [min_indel_cnt / 2, max_indel_cnt * 2],
          'y-ticks': [min_indel_cnt, max_indel_cnt],
          'y-t-labels': [min_indel_cnt, max_indel_cnt]}
  }
  return axes_specs


def decorate_axes(ax, axes_spec):
  """Decorate axes as needed"""
  # Set axes x-lims. We operate by columns here
  fxl = axes_spec['x-lim']
  x_lim = {'DEL': [fxl[0], -1], 'REF/SNP': [-1, 1], 'INS': [1, fxl[1]]}
  fxt = axes_spec['x-ticks']
  x_ticks = {'DEL': [t for t in fxt if t < 0], 'REF/SNP': [-0.5, 0.5], 'INS': [t for t in fxt if t > 0]}
  x_tick_labels = {'DEL': [abs(x) for x in x_ticks['DEL']], 'REF/SNP': ['Ref', 'SNP'], 'INS': x_ticks['INS']}

  for k1 in ['DEL', 'REF/SNP', 'INS']:
    for k2 in ['A', 'Ad', 'B', 'C']:
      if k2 not in axes_spec: continue
      y_lim, y_ticks = axes_spec[k2]['y-lim'], axes_spec[k2]['y-ticks']
      y_tick_labels = axes_spec[k2].get('y-t-labels', y_ticks)
      ak = k2 + '-' + k1
      if k2 not in ['A', 'Ad']:
        ax[ak].set_yscale('symlog', nonposy='clip', subsy=[])
      plt.setp(ax[ak],
               xlim=x_lim[k1], xticks=x_ticks[k1], xticklabels=[] if k2 != 'C' else x_tick_labels[k1],
               ylim=y_lim,
               yticks=y_ticks,
               yticklabels=[] if k1 != 'DEL' else y_tick_labels)
      if ak == 'C-REF/SNP':
        ax[ak].set_xticklabels(x_tick_labels[k1], rotation=90)

      ax[ak].xaxis.grid('true')  # The vertical guide lines

      if k2 in ['A', 'Ad']:
        ax[ak].yaxis.grid('true')  # The horizontal guide lines

      # Ticks and axes
      for p in ['left', 'top', 'right']:
        if k1 == 'DEL' and p == 'left': continue
        ax[ak].spines[p].set_visible(False)

      ax[ak].minorticks_off()
      ax[ak].tick_params(bottom='off', left='off', top='off', right='off')
      if k2 == 'C':
        ax[ak].tick_params(direction='out', length=3, bottom='on', pad=3)
      if k1 == 'DEL':
        ax[ak].tick_params(direction='out', length=3, left='on', pad=2)

  ax['A-REF/SNP'].title.set(text=axes_spec['title'])
  ax['A-DEL'].legend(title='|d|', loc='lower left', fontsize=8)
  if 'Ad-DEL' in ax: ax['Ad-DEL'].legend(loc='lower left', fontsize=8)

  # X-labels
  ax['C-DEL'].set_xlabel('Deletions (bp)')
  ax['C-INS'].set_xlabel('Insertions (bp)')

  # Y-labels
  ax['A-DEL'].set_ylabel('% reads\ncorrectly aligned')
  if 'Ad-DEL' in ax: ax['Ad-DEL'].set_ylabel('Difference')
  ax['B-DEL'].set_ylabel('Read\ncount')
  ax['C-DEL'].set_ylabel('Variation\ncount')


def plot_alignment_accuracy(read_counts, d_error_cat, color='k', lw=2, ax={}):
  """axs is a dict of axes A-DEL, A-REF/SNP and A-INS"""
  #pc_raw = cat_read_counts['correct'] / cat_read_counts['total'].astype(float) * 100
  #pc = ss.medfilt(np.ma.masked_invalid(pc_raw), kernel_size=kernel_size) if kernel_size > 2 else pc_raw

  rc = read_counts

  v_pc = rc['v_r_cnt'][:, :d_error_cat[1] + 1].sum(axis=1) / rc['v_r_cnt'].sum(axis=1).astype(float) * 100

  # Indels
  idx = np.where(rc['v_len'] < 0)
  ax['A-DEL'].plot(rc['v_len'][idx], v_pc[idx], color=color, lw=lw, label=d_error_cat[0], alpha=0.71)
  idx = np.where(rc['v_len'] > 0)
  ax['A-INS'].plot(rc['v_len'][idx], v_pc[idx], color=color, lw=lw, label=d_error_cat[0], alpha=0.71)

  # SNP
  idx = np.where(rc['v_len'] == 0)[0]
  x, y = 0.5, v_pc[idx]
  ax['A-REF/SNP'].plot(x, y, color=color, marker='o', ms=2, mec=color, alpha=0.71)
  ax['A-REF/SNP'].plot(x, y, color=color, marker='_', lw=1.5, alpha=0.71)

  # ref reads
  ref_pc = rc['ref_r_cnt'][:d_error_cat[1] + 1].sum() / float(rc['ref_r_cnt'].sum()) * 100
  ax['A-REF/SNP'].plot(-0.5, ref_pc, color=color, marker='o', ms=2, mec=color, alpha=0.71)
  ax['A-REF/SNP'].plot(-0.5, ref_pc, color=color, marker='_', lw=1.5, alpha=0.71)


def plot_read_count(read_counts, color='k', lw=2, ax={}):
  rc = read_counts

  indel_rc = rc['v_r_cnt'].sum(axis=1)

  # Indels
  idx = np.where(rc['v_len'] < 0)
  ax['B-DEL'].plot(rc['v_len'][idx], indel_rc[idx], color=color, lw=lw, alpha=0.71)
  idx = np.where(rc['v_len'] > 0)
  ax['B-INS'].plot(rc['v_len'][idx], indel_rc[idx], color=color, lw=lw, alpha=0.71)

  # SNP
  idx = np.where(rc['v_len'] == 0)[0]
  x, y = 0.5, indel_rc[idx]
  ax['B-REF/SNP'].plot(x, y, color=color, marker='o', ms=2, mec=color, alpha=0.71)
  ax['B-REF/SNP'].plot(x, y, color=color, marker='_', lw=1.5, alpha=0.71)

  # ref reads
  ref_rc = rc['ref_r_cnt'].sum()
  ax['B-REF/SNP'].plot(-0.5, ref_rc, color=color, marker='o', ms=2, mec=color, alpha=0.71)
  ax['B-REF/SNP'].plot(-0.5, ref_rc, color=color, marker='_', lw=1.5, alpha=0.71)


def plot_var_count(read_counts, color='k', lw=2, ax={}):
  rc = read_counts

  indel_vc = rc['v_cnt']

  # Indels
  idx = np.where(rc['v_len'] < 0)
  ax['C-DEL'].plot(rc['v_len'][idx], indel_vc[idx], color=color, lw=lw, alpha=0.71)
  idx = np.where(rc['v_len'] > 0)
  ax['C-INS'].plot(rc['v_len'][idx], indel_vc[idx], color=color, lw=lw, alpha=0.71)

  # SNP
  idx = np.where(rc['v_len'] == 0)[0]
  x, y = 0.5, indel_vc[idx]
  ax['C-REF/SNP'].plot(x, y, color=color, marker='o', ms=2, mec=color, alpha=0.71)
  ax['C-REF/SNP'].plot(x, y, color=color, marker='_', lw=1.5, alpha=0.71)


if __name__ == '__main__':
  cli()