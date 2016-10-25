"""Example read model

TODO: Figure out the logistics of having diffrent models and loading GC-bias, sequencing error etc.

"""
import pickle

import numpy as np


SEED_MAX = (1 << 32) - 1  # Used for seeding rng


def read_model_params(model, diploid_coverage=30.0):
  """

  :param model
  :param diploid_coverage: Coverage we expect from a diploid genome
  :return: dict of params,
           passes -> how many times we should call generate reads for each region
  """
  rlen = model['mean_rlen']
  # coverage = read count * read len / genome len
  # coverage = P * read len P = expected number of reads per base
  # P = p * passes    p = probability of read per pass
  # One template for two reads
  # p = P / (2 * passes)
  # p = coverage / (2 * rlen * passes)
  p = 1.0
  passes = 1
  while p > 0.1:
    passes *= 2
    p = 0.5 * diploid_coverage / (2 * rlen * passes)

  return {
    'diploid_coverage': diploid_coverage,
    'p': p,
    'passes': passes,
    'rlen': rlen,
    'cum_tlen': model['cum_tlen'],
    'cum_bq_mat': model['cum_bq_mat']
  }


def generate_reads(model, p_min, p_max, seed=7):
  """

  :param model: dict as returned by read_model_params
  :param p_min: Start of reads
  :param p_max: End of reads
  :param region: ('1', 10001, 103863906) Like the BED file
  :param seed:
  :return:
  """
  if not (0 <= seed <= SEED_MAX):
    raise ValueError('Seed value {} is out of range 0 - {}'.format(seed, SEED_MAX))

  seed_rng = np.random.RandomState(seed)
  tloc_rng, tlen_rng, shuffle_rng, file_order_rng = \
    [np.random.RandomState(s) for s in seed_rng.randint(SEED_MAX, size=4)]

  return _reads_for_template_in_region(
        model, file_order_rng,
        *_templates_for_region(
          p_min, p_max, model, tloc_rng, tlen_rng, shuffle_rng))


def _templates_for_region(p_min, p_max, model, tloc_rng, tlen_rng, shuffle_rng):
  # TODO: GC-bias goes here
  p, rlen = model['p'], model['rlen']
  est_block_size = int((p_max - p_min) * p * 1.2)
  ts = tloc_rng.geometric(p=p, size=est_block_size).cumsum() + p_min + 1
  shuffle_rng.shuffle(ts)
  tl = np.searchsorted(model['cum_tlen'], tlen_rng.rand(ts.shape[0]))
  # tl = (tlen_rng.randn(ts.shape[0]) * tlen_std + tlen).astype(int)
  tl.clip(rlen)
  te = ts + tl
  idx = (te < p_max)
  return ts[idx], te[idx]


def _reads_for_template_in_region(model, file_order_rng, ts, te):
  """

  :param rlen: Fixed - Illumina reads
  :param file_order_rng: Decides which one of the pair (F or R) comes first in the file (or goes in file1 of the pair)
  :param ts: Start of template
  :param te: End of template
  :return:
  """
  rlen = model['rlen']

  r0l = np.full(ts.size, rlen, dtype=np.uint32)
  r1l = np.full(ts.size, rlen, dtype=np.uint32)

  r0fo = file_order_rng.randint(2, size=ts.size, dtype='i1')
  r1fo = 1 - r0fo

  r0p = ts
  r1p = te - rlen

  return [
    {
      'file_order': r0fo,
      'pos': r0p,
      'len': r0l
    },
    {
      'file_order': r1fo,
      'pos': r1p,
      'len': r1l
    }
  ]


def describe_model(model, figfile):
  """Plot a few panels describing what the model looks like

  :param model:
  :return:
  """
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt

  fig = plt.figure(figsize=(6, 15))
  plt.subplots_adjust(bottom=0.05, top=0.97)

  ax = plt.subplot(4,1,1)
  ax.text(0.01, 0.99, 'Model description:\n\n' + model['model_description'], va='top', wrap='True')
  plt.setp(ax, xlim=(0, 1), ylim=(0, 1), xticks=[], yticks=[])

  ax = plt.subplot(4,1,2)
  plot_template_length_distribution(model, ax, plt)

  ax = plt.subplot(4,1,3)
  plot_BQ_heatmap(model, 0, ax, plt)

  ax = plt.subplot(4,1,4)
  plot_BQ_heatmap(model, 1, ax, plt)

  plt.savefig(figfile)


def plot_template_length_distribution(model, ax, plt):
  ax.plot(model['tlen'])
  ax.text(len(model['tlen']), plt.getp(ax, 'ylim')[1], 'n={}'.format(model['r_cnt']), ha='right', va='top')
  plt.setp(ax, xlabel='Template length (bp)', ylabel='Read count')


def plot_BQ_heatmap(model, mate, ax, plt):
  plt.imshow(model['bq_mat'][mate, :65, :model['max_rlen']].T,
             extent=(1, model['max_rlen'], 1, 65),
             aspect='auto',
             origin='lower', cmap=plt.cm.gray_r)
  ax.text(plt.getp(ax, 'xlim')[1], plt.getp(ax, 'ylim')[1], 'Mate {}'.format(mate + 1),
          ha='right', va='top')
  plt.setp(ax, xlabel='Position on read (bp)', ylabel='BQ value')