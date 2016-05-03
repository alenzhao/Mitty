import itertools
import cPickle

import click
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.signal as ss

import mitty.benchmarking.creed as creed

import logging
logger = logging.getLogger(__name__)


variant_categories = [
  ('200 < DEL', float('-inf'), -200),  # (key, i, j)  =>  i <= L < j
  ('100 < DEL <= 200', -200, -100),
  ('50 < DEL <= 100', -100, -50),
  ('20 < DEL <= 50', -50, -20),
  ('0 < DEL <= 20', -20, 0),
  ('SNP', 0, 1),
  ('0 < INS <= 20', 1, 21),
  ('20 < INS <= 50', 21, 51),
  ('50 < INS <= 100', 51, 101),
  ('100 < INS <= 200', 101, 201),
  ('200 < INS', 201, float('+inf')),
]


@click.command()
@click.version_option()
@click.argument('indelpkl', type=click.Path(exists=True))
@click.argument('outsuffix', type=click.Path())
@click.option('--title', help='Title', default='Aligner accuracy')
def cli(indelpkl, outsuffix, title):
  """pickled data -> plot:

  \b
  Given a pkl file with alignment accuracy parametrized by the sample indel lengths,
  generate a plot of alignment accuracy with |d|."""
  read_counts = cPickle.load(open(indelpkl, 'r'))

  for sl in ['known', 'novel']:
    plot_fig(creed.slice_read_counts(read_counts, sl), '{}-{}'.format(title, sl))
    out_name = '{}-{}'.format(sl, outsuffix)
    plt.savefig(out_name)
    print('Saved {} variant plot in {}'.format(sl, out_name))


# TODO: fix indexes for variant ranges. Plot is OK for first pass
def plot_fig(read_counts, title):

  indel_range = (read_counts['v_r_cnt'].shape[0] - 1 )/2
  available_variant_categories = [
    vc for vc in variant_categories
    if (0 <= vc[1] + indel_range < 2 * indel_range + 1) or (0 < vc[2] + indel_range <= 2 * indel_range + 1)
  ]

  # We share the x-axis
  f, axarr = plt.subplots(len(available_variant_categories) + 1, sharex=True, sharey=True, figsize=(5,20))  # +1 because of Ref plot
  f.subplots_adjust(hspace=0.02, bottom=0.05, top=0.98, right=0.87)
  d_range = read_counts['ref_r_cnt'].shape[0]

  axarr[0].plot(100.0 * read_counts['ref_r_cnt'].cumsum() / read_counts['ref_r_cnt'].sum())
  plt.setp(axarr[0], xlim=[-5, 100], ylabel='% Correct',ylim=[0, 110], yticks=[50, 60, 70, 80, 90, 100])
  axarr[0].text(0, 10, 'REF', backgroundcolor=(.7,.7,.7))

  for n, spec in enumerate(available_variant_categories):
    i = max(0, spec[1] + indel_range)
    j = max(0, min(spec[2] + indel_range, 2 * indel_range + 1))
    # if j < 1: continue
    axarr[n + 1].plot(100.0 * read_counts['v_r_cnt'][i:j, :].sum(axis=0).cumsum() / read_counts['v_r_cnt'][i:j, :].sum())
    axarr[n + 1].text(0, 10, spec[0], backgroundcolor=(.7,.7,.7))

  for n in range(len(available_variant_categories) + 1):
    axarr[n].yaxis.tick_right()
    # axarr[n].yaxis.set_label_position('right')
    axarr[n].yaxis.grid(True)
    axarr[n].xaxis.grid(True)
    axarr[n].tick_params(axis='y', labelsize=8)

  plt.setp(axarr[n], xlabel='|d|')