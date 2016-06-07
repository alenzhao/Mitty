"""In the final form, this will collate all the benchmarking data (aligner + VC) and produce summary reports.

V0 will

a) Collate summary VC data (from the eval.csv) and create a pandas panel from them.
b) Create an html webpage Leaderboard that arranges aligners by VC performance.
   There are multiple leader boards generated based on the different metadata criteria
   In the current case they are:

   corruption (x2), variant caller (x3) and graph (x4)

Input is a csv file with every row as follows

  run
  tool
  sample
  graph
  corrupt
  read_len
  template_len
  eval.csv file path
  aligner_analysis task url
  variant caller task url
"""
import csv
import StringIO
import base64

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import click
import numpy as np
import pandas as pd

import logging
logger = logging.getLogger(__name__)


__metadata_keys__ = [
  'aligner',
  'caller',
  'run',
  'sample',
  'graph',
  'corrupt',
  'read_len',
  'template_len',
]


@click.group()
@click.version_option()
def cli():
  pass


@cli.command('variant-leaderboard')
@click.argument('incsv', type=click.Path(exists=True))
@click.argument('outcsv', type=click.Path())
@click.argument('outhtml', type=click.Path())
def vc_leaderboard(incsv, outcsv, outhtml):
  """csv -> leaderboard + pandas Panel in hdf5

  \b

  Create an html webpage Leaderboard that arranges aligners by VC performance.

  input csv file columns are

  run
  aligner
  sample
  graph
  corrupt
  read_len
  template_len
  caller
  eval_csv_path
  aligner_analysis task url
  variant caller task url

  The output data frame is saved as a CSV and contains this data combined together.

  See http://cruncher-new:8888/notebooks/kghose/LeaderBoard/data_analysis_example.ipynb#
  for how to analyse the data

  """
  in_data = [{k.strip(): v.strip() for k, v in row.items()} for row in csv.DictReader(open(incsv, 'r'))]
  fname_list = [row['eval_csv_path'] for row in in_data]
  meta_row = in_data
  pan = combine_evaluation_csvs(fname_list, meta_row)
  pan.to_csv(outcsv)
  # Read back with:
  # pan = pd.read_csv('aggregated-vc-data-2016-06-03T13-41-16.227566.csv', index_col=[0, 1, 2, 3, 4, 5, 6, 7], skipinitialspace=True, header=[0,1])

  with open(outhtml, 'w') as fp:
    fp.write(html_leaderboard(pan))


def combine_evaluation_csvs(fname_list, meta_row):
  data = [pr_by_variant_type(load_evaluation_csv(fname)) for fname in fname_list]
  index_tuples = [tuple([meta[key] for key in __metadata_keys__]) for meta in meta_row]
  return pd.DataFrame(data, index=pd.MultiIndex.from_tuples(index_tuples, names=__metadata_keys__))


def load_evaluation_csv(fname):
  use_cols = ['Type', 'Subtype', 'Genotype', 'Filter',
              'METRIC.Recall', 'METRIC.Precision', 'METRIC.F1_Score',
              'TRUTH.TOTAL', 'TRUTH.TP', 'TRUTH.FN',
              'QUERY.TOTAL', 'QUERY.TP', 'QUERY.FP']
  index_cols = ['Type', 'Subtype', 'Genotype', 'Filter']
  df = pd.read_csv(fname, usecols=use_cols, index_col=index_cols)
  return df


def pr_by_variant_type(df):
  """Extract the P/R for the built in variant types"""
  indel_size = [('INDEL', 'D16_PLUS'), ('INDEL', 'D6_15'), ('INDEL', 'D1_5'),
                ('SNP', '*'),
                ('INDEL', 'I1_5'), ('INDEL', 'I6_15'), ('INDEL', 'I16_PLUS')]

  pr = np.concatenate([df.xs(isize + ('*', 'ALL'), level=('Type', 'Subtype', 'Genotype', 'Filter')).loc[:, ('METRIC.Precision', 'METRIC.Recall')].values for isize in indel_size])
  row_names = [i[1] if i[1] != '*' else 'SNP' for i in indel_size]
  new_df = pd.DataFrame(pr, columns=['Precision', 'Recall'], index=row_names)
  new_df['F'] = 2 * new_df['Precision'] * new_df['Recall'] / (new_df['Precision'] + new_df['Recall'])

  return new_df.T.unstack()


html_header = """
<!DOCTYPE html>
<html>
<head>
<style>
  body {
   font-family: Gill Sans, sans-serif;
   font-size: 9pt;
   padding-left: 20px;
   padding-right: 20px;
  }

  table {
      border-collapse: collapse;
      width: 100%;
  }

  th, td {
      padding: 8px;
      text-align: left;
      border-bottom: 1px solid #ddd;
      border-left: 1px solid #ddd;
  }
</style>
</head>
<body>
"""

html_footer = """</body>"""


def html_leaderboard(pan):
  """This is a set of plots organized as a nested table"""
  html = []
  for caller in ['SBG-JB', 'GATK-3.5', 'FB']:
    html += create_table_for_one_caller(pan, caller)

  return html_header + '\n'.join(html) + html_footer


def create_table_for_one_caller(pan, caller='SBG-JB'):
  graphs = ['None', 'S', 'G0', 'G1']
  variant_classes = ['SNP', 'I1_5', 'D1_5', 'I6_15', 'D6_15', 'I16_PLUS', 'D16_PLUS']

  html = ['<h2>Caller: {}</h2>'.format(caller)]
  html += ['<table>']
  html += ['<tr>'] + ['<th></th>'] + ['<th colspan=2>graph={}</th>'.format(g) for g in graphs] + ['</tr>']
  html += ['<tr>'] + ['<th>Variant<br>Class</th>'] + ['<th>corrupt={}</th>'.format(c) for g in graphs for c in ['No', 'Yes']] + ['</tr>']
  for variant_class in variant_classes:
    html += ['<tr>'] + ['<th>{}</th>'.format(variant_class)]
    for graph in graphs:
      for corrupt in ['No', 'Yes']:
        html += ['<td>{}</td>'.format(draw_pr_plot(pan, caller, graph, corrupt, variant_class))]
  html += ['</table>']
  return html


def draw_pr_plot(pan, caller, graph, corrupt, variant_class):
  this_slice = {
    'key': (1, 'S', graph, corrupt, 100, 500, caller),
    'level': ('run', 'sample', 'graph', 'corrupt', 'read_len', 'template_len', 'caller')}

  fig = plt.figure(figsize=(4,8))
  fig.subplots_adjust(bottom=0.25, right=0.99, top=0.90)
  ax = fig.add_subplot(111)

  pan.xs(**this_slice)[variant_class].sort_index(ascending=False).ix[:, :'Recall'].plot(ax=ax, rot=90, lw=4, style='o-')
  plt.title('Graph: {}\nVariant: {}\nCorrupt: {}'.format(graph, variant_class, corrupt))

  imgdata = StringIO.StringIO()
  plt.savefig(imgdata, format='png')
  imgdata.seek(0)  # rewind the data
  plt.close()

  return '<img src="data:image/png;base64,{}" width="200px"></img>'.format(base64.b64encode(imgdata.buf))


if __name__ == '__main__':
  cli()

