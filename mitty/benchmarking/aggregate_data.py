import cPickle
import csv

import click
import plotly.offline as py
import plotly.graph_objs as go

import mitty.benchmarking.creed as creed


import logging
logger = logging.getLogger(__name__)


@click.command()
@click.version_option()
@click.argument('incsv', type=click.Path(exists=True))
@click.argument('outcsv', type=click.Path())
@click.argument('outhtml', type=click.Path())
@click.option('--indel-range', default=100, help='Range of indels to plot')
@click.option('--median-filter-size', default=0, help='Size of median filter window to smooth plots')
@click.option('--title', help='Title', default='Aligner accuracy')
def cli(incsv, outcsv, outhtml, indel_range, median_filter_size, title):
  """csv -> indel-pkls -> plot
                       -> csv file with data points (can be used for further plotting)
  \b
  Given a csv file with a list of indel-pkl files + metadata generate an interactive plot

  input csv file columns are

  tool-version-string
  sample
  graph
  read-corrupted
  read-len
  indel-pkl file name

  The output csv contains the same initial columns followed by data columns with the indel size (or Ref)
  marked out in the header.

  -500, -499, -498 .... Ref, SNP, 1, 2, 3, 4 ..... 500
  """
  median_filter_size += (median_filter_size + 1) % 2  # Ensure median_filter_size is odd

  column_names = column_names_from_indel_range(indel_range)

  with csv.DictWriter(open(outcsv, 'w'), fieldnames=column_names) as outcsv_fp:
    for row in csv.DictReader(open(incsv, 'r')):
      out_row = {k: row[k] for k in ['tool_version', 'sample', 'graph', 'reads_corrupt', 'read_len', 'reads_paired']}
      out_row = dict(out_row, )




  read_counts = cPickle.load(open(indelpkl, 'r'))

  for sl in ['known', 'novel']:
    plot_fig(creed.slice_read_counts(read_counts, sl), indel_range, median_filter_size, '{}-{}'.format(title, sl))
    out_name = '{}-{}'.format(sl, outsuffix)
    plt.savefig(out_name)
    print('Saved {} variant plot in {}'.format(sl, out_name))


def column_names_from_indel_range(indel_range):
  return ['tool_version', 'sample', 'graph', 'reads_corrupt', 'read_len', 'reads_paired'] \
         + [str(x) for x in range(-indel_range, 0)] + ['Ref', 'SNP'] \
         + [str(x) for x in range(1, indel_range + 1)]
