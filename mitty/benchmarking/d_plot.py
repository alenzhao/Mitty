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


# variant_categories = [
#   ('200 < DEL', float('-inf'), -200),  # (key, i, j)  =>  i <= L < j
#   ('100 < DEL <= 200', -200, -100),
#   ('50 < DEL <= 100', -100, -50),
#   ('20 < DEL <= 50', -50, -20),
#   ('0 < DEL <= 20', -20, 0),
#   ('SNP', 0, 1),
#   ('0 < INS <= 20', 1, 21),
#   ('20 < INS <= 50', 21, 51),
#   ('50 < INS <= 100', 51, 101),
#   ('100 < INS <= 200', 101, 201),
#   ('200 < INS', 201, float('+inf')),
# ]


variant_categories = [
  ('100 < DEL', float('-inf'), -100),  # (key, i, j)  =>  i <= L < j
  ('50 < DEL <= 100', -100, -50),
  ('20 < DEL <= 50', -50, -20),
  ('0 < DEL <= 20', -20, 0),
  ('SNP', 0, 1),
  ('0 < INS <= 20', 1, 21),
  ('20 < INS <= 50', 21, 51),
  ('50 < INS <= 100', 51, 101),
  ('100 < INS', 101, float('+inf')),
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


def plot_fig(read_counts, title):
  indel_range = (read_counts['v_r_cnt'].shape[0] - 1 )/2
  d = range(-creed.MAX_D_ERROR, creed.MAX_D_ERROR + 1)

  plt.figure(figsize=(10, 5))
  y_xneg = 100.0 * read_counts['v_r_cnt'][:, creed.MAX_D_ERROR::-1].cumsum(axis=1)[:, ::-1]
  y_xpos = 100.0 * read_counts['v_r_cnt'][:, creed.MAX_D_ERROR:].cumsum(axis=1)  # Need to duplicate d = 0 point for cumulative sum
  y = np.concatenate((y_xneg, y_xpos[:, 1:-1]), axis=1)  # get rid of duplicate d = 0 point, get rid of d = MAX_D + 1 which is unmapped and reads on wrong chrom
  Z = y/np.clip((y[:, 0] + y[:, -1] - y[:, creed.MAX_D_ERROR])[:, None], 1, 1e10)
  idx = np.nonzero(Z[:, creed.MAX_D_ERROR])[0]
  clim = (Z[idx, creed.MAX_D_ERROR].min(), Z.max()) if idx.size else (0, 1.0)
  plt.imshow(Z, extent=(d[0], d[-1], -indel_range, indel_range),
              cmap=plt.cm.gray, clim=clim,
              interpolation='none',
              origin='lower', aspect='auto')
  plt.colorbar()

  plt.setp(plt.gca(), xlim=[-100, 100])
  plt.xlabel('d')
  plt.ylabel('INDEL size')
  plt.subplots_adjust(left=0.1, right=0.98)
  plt.suptitle(title)