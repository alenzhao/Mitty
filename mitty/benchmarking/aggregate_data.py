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
  There are multiple leader boards generated based on the different metadata criteria
  In the current case they are:

    corruption (x2), variant caller (x3) and graph (x4)

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

  The output data frame contains this data combined together.
  """
  in_data = [{k.strip(): v.strip() for k, v in row.items()} for row in csv.DictReader(open(incsv, 'r'))]
  fname_list = [row['eval_csv_path'] for row in in_data]
  meta_row = in_data
  pan = combine_evaluation_csvs(fname_list, meta_row)
  pan.to_csv(outcsv)

  with open(outhtml, 'w') as fp:
    fp.write(nicely_formatted_leader_board(leader_board(pan).to_html()))



def combine_evaluation_csvs(fname_list, meta_row):
  data = [pr_by_variant_type(load_evaluation_csv(fname)) for fname in fname_list]
  index_tuples = [tuple([meta[key] for key in __metadata_keys__]) for meta in meta_row]
  #return pd.concat(data, axis=0, keys=index_tuples, names=__metadata_keys__)
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


def default_scoring_function(df):
  """This simply takes the mean F number across all indel categories and ranks the runs accordingly."""
  return df.xs('F', level=1, axis=1).mean(axis=1).sort_values(ascending=False)


def leader_board(df, scoring_fun=None, criteria=None):
  """Given a data frame, a scoring function and a set of criteria, organize the data frame to show us the aligner names
  and metrics sorted according to the scoring function


  :param df:  data frame
  :param scoring_fun: function that takes in a data frame and returns a ranking of rows
                      see the default scoring function for an example
  :param criteria:
  :return:
  """
  #return df.reindex(df.xs('F', level=1, axis=1).mean(axis=1).sort_values(ascending=False).index)
  scoring_fun = scoring_fun or default_scoring_function
  return df.reindex(scoring_fun(df).index).T


def nicely_formatted_leader_board(html_table):
  """Given the raw table from Pandas combine it with a style sheet. Remember to set the class to be leaderboard"""
  top = """<!DOCTYPE html>
<html>
<head>
<style>
  body {
   font-family: Gill Sans, sans-serif;
   padding-left: 100px;
   padding-right: 100px;
  }

  table {
      border-collapse: collapse;
      width: 100%;
  }

  th, td {
      padding: 8px;
      text-align: left;
      border-bottom: 1px solid #ddd;
  }
</style>
</head>
<body>

<h2>Aligner/Caller leaderboard</h2>
"""
  bottom = """
<b>The leaderboard rank is based on the average F value across all variant classes</b>
</body>
</html>"""
  return top + html_table + bottom


def plot_me():
  f = 2 * pr[:,0] * pr[:,1]/(pr[:,0] + pr[:,1])
  plt.plot(pr[:,0], pr[:,1], 'b')
  plt.plot(pr[:3,0], pr[:3,1], 'ko', ms=10)
  plt.plot(pr[3,0], pr[3,1], 'ks', ms=10)
  plt.plot(pr[4:,0], pr[4:,1], 'k+', ms=10)
  plt.axis('scaled')
  plt.show()





def pr_computation(df):
  """Given the df from a variant calling run, compute the P/R for the different categories available."""
  x = df  # .xs('SNP', level='Type')
  p = x['TRUTH.TP'] / (x['TRUTH.TP'] + x['QUERY.FP'])
  r = x['TRUTH.TP'] / (x['TRUTH.TP'] + x['TRUTH.FN'])
  f = 2 * p * r / (p + r)


if __name__ == '__main__':
  cli()

