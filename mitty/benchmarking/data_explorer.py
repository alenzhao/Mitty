"""
The following 2D slices into this multi-dimensional beast are available. y-axis is always P+R

x-axis is over

aligner
caller
variant_category
graph
corrupt



"""
import pandas as pd
import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec


def load_data(fname):
  #pan = pd.read_csv('/Users/kghose/TestData/variant-calling/aggregated-vc-data-2016-06-03T13-41-16.227566.csv', index_col=[0, 1, 2, 3, 4, 5, 6, 7], skipinitialspace=True, header=[0,1])
  pan = pd.read_csv(fname, index_col=[0, 1, 2, 3, 4, 5, 6, 7], skipinitialspace=True, header=[0,1])
  df = pan.stack(0)[['Precision', 'Recall']]
  i2 = df.index
  names = list(i2.names)
  names[-1] = 'variant_category'
  i2.names = names
  df.index = i2
  reordered_names = ['run', 'sample', 'read_len', 'template_len', 'corrupt', 'aligner', 'caller', 'graph', 'variant_category']
  df2 = df.reorder_levels(reordered_names)
  return df2


def setup_axes():
  fig = plt.figure(figsize=(12, 8))
  fig.subplots_adjust(bottom=0.05, left=0.1, right=0.99, top=0.98, hspace=0.3, wspace=0.05)

  ax = {}
  ax['aligner'] = plt.subplot2grid((2, 5), (0, 0), colspan=3)
  ax['caller'] = plt.subplot2grid((2, 5), (0, 3), colspan=2, sharey=ax['aligner'])
  ax['variant_category'] = plt.subplot2grid((2, 5), (1, 0), colspan=2)
  ax['graph'] = plt.subplot2grid((2, 5), (1, 2), colspan=2, sharey=ax['variant_category'])
  ax['corrupt'] = plt.subplot2grid((2, 5), (1, 4), sharey=ax['variant_category'])
  return fig, ax


def setup_plot_args():
  return {
    'aligner': {'rot': 90, 'lw': 2, 'style': 'o-'},
    'caller': {'rot': 90, 'lw': 2, 'style': 'o-'},
    'variant_category': {'lw': 2, 'style': 'o-'},
    'graph': {'lw': 2, 'style': 'o-'},
    'corrupt': {'lw': 2, 'style': 'o-'}
  }


def setup_index_order(pr_data):
  return {
    'aligner': pr_data.index.levels[pr_data.index.names.index('aligner')].sort_values(ascending=False),
    'caller': pr_data.index.levels[pr_data.index.names.index('caller')].sort_values(ascending=False),
    'variant_category': pd.Index(['D16_PLUS', 'D6_15', 'D1_5', 'SNP', 'I1_5', 'I6_15', 'I16_PLUS']),
    'graph': pd.Index(['None', 'S', 'G0', 'G1']),
    'corrupt': pd.Index(['No', 'Yes'])
  }


def get_xp_slice(df, slice_dict, this_slice_key, s_index):
  """

  :param df:
  :param slice_dict: {k: v} pairs indicating the center of the slice e.g.
                     {'run': 1, 'sample': 'S', .... }
  :param this_slice_key: key of the dimension we are expanding here
  :param s_index: the index sorted the way we want
  :return: a data frame with just this slice
  """
  keys = slice_dict.keys()
  this_slice = {
    'key': tuple([slice_dict[k] for k in keys if k != this_slice_key]),
    'level': tuple([k for k in keys if k != this_slice_key])
  }
  return df.xs(**this_slice).reindex(s_index).ix[:, :'Recall']


def xy_slice(df, slice_dict, this_slice_key, s_index, ax, **plot_kwargs):
  """Plot the data frame using the given slice

  :param df:
  :param slice_dict: {k: v} pairs indicating the center of the slice e.g.
                     {'run': 1, 'sample': 'S', .... }
  :param this_slice_key: key of the dimension we are expanding here
  :param s_index: the index sorted the way we want
  :param ax: the axes we want to plot on
  :param plot_kwargs: anything to be passed to plot, e.g. rot=90, lw=4, style='o-'
  :return: No return value - just a plot

  """
  keys = slice_dict.keys()
  this_slice = {
    'key': tuple([slice_dict[k] for k in keys if k != this_slice_key]),
    'level': tuple([k for k in keys if k != this_slice_key])
  }
  sl = df.xs(**this_slice).reindex(s_index).ix[:, :'Recall']
  sl.plot(ax=ax, fontsize=9, xlim=(-0.5, len(s_index) - 0.5), xticks=range(len(s_index)), **plot_kwargs)
  #df.xs(**this_slice).reindex(s_index).ix[:, :'Recall'].plot(ax=ax, fontsize=9, **plot_kwargs)
  ax.axvline(s_index.get_loc(slice_dict[this_slice_key]), color='r', linestyle='--', lw=2)


def plot_panel(pr_data, sl_d, ax, io, pl_args):
  for k in ['aligner', 'caller', 'variant_category', 'graph', 'corrupt']:
    xy_slice(pr_data, sl_d, k, io[k], ax[k], **pl_args[k])


class DataExplorer:
  def __init__(self, fname):
    self.pr_data = load_data(fname)

    # The starting slice
    self.sl_d = {
      'run': 1,
      'sample': 'S',
      'read_len': 100,
      'template_len': 500,
      'corrupt': 'Yes',
      'aligner': '0.8.8',
      'caller': 'SBG-JB',
      'graph': 'None',
      'variant_category': 'SNP'
    }

    self.fig, self.ax = setup_axes()
    self.io = setup_index_order(self.pr_data)
    self.pl_args = setup_plot_args()

    self.replot()

    cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)

  def replot(self):
    for ax in self.ax.values():
      ax.clear()
    plot_panel(self.pr_data, self.sl_d, self.ax, self.io, self.pl_args)
    plt.draw()

  def find_axes(self, event):
    return next((k for k in self.ax.keys() if self.ax[k] == event.inaxes), None)

  def onclick(self, event):
    clicked_ax = self.find_axes(event)
    if clicked_ax is None:
      return
    self.sl_d[clicked_ax] = self.io[clicked_ax].values[int(round(event.xdata))]

    self.replot()
    # print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
    #       (event.button, event.x, event.y, event.xdata, event.ydata))

#de = DataExplorer('hi')