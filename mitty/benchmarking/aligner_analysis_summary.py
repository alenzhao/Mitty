import json
from collections import OrderedDict
import logging
import base64

import click
import numpy as np

logger = logging.getLogger(__name__)


__html_template__ = \
"""
<html>
<head>
  <title>{}</title>
  <style>
    body {{
     font-family: Gill Sans, sans-serif;
     padding-left: 100px;
     padding-right: 100px;
    }}
    .metadata {{
      border-collapse: collapse;
      border: 1px solid black;
    }}
    .metadata tr:nth-child(even) {{background-color: #dfdfdf;}}
    .metadata  tr:nth-child(odd) {{background-color: #cfcfcf;}}
    .metadata th, td {{
      padding-right: 15px;
      padding: 5px;
      text-align: left;
      border: 1px solid black;
    }}


    .accuracy {{
      border-collapse: collapse;
      border: 1px solid black;
      font-size: 12px;
      font-weight: normal;
    }}
    .accuracy th {{
        border: 1px solid black;
    }}
    .accuracy td {{
        padding-right: 15px;
        padding: 5px;
        text-align: right;
        border: 1px solid black;
    }}

    div {{
      background-color: yellow;
      padding: 2px;
      display: none;
    }}
    span:hover {{
      background-color: red;
    }}
    span:hover + div {{
      display: table;
      position: absolute;
    }}

  </style>
</head>
<body>
{}
</body>
</html>
"""


@click.command()
@click.argument('tool')
@click.option('--sample', default='Ref')
@click.option('--graph', default='-')
@click.option('--db-summary', type=click.Path(exists=True))
@click.option('--known-d-plot', type=click.Path(exists=True))
@click.option('--novel-d-plot', type=click.Path(exists=True))
@click.argument('mismatcsv', type=click.Path(exists=True))
@click.argument('summaryjson', type=click.Path(exists=True))
@click.argument('novelplot', type=click.Path(exists=True))
@click.argument('knownplot', type=click.Path(exists=True))
@click.argument('mqplot', type=click.Path(exists=True))
@click.argument('circleplot', type=click.Path(exists=True))
@click.argument('matrixplot', type=click.Path(exists=True))
@click.argument('htmlout', type=click.Path())
def cli(tool, sample, graph, mismatcsv, summaryjson, novelplot, knownplot, mqplot, circleplot, matrixplot,
        db_summary, known_d_plot, novel_d_plot,
        htmlout):
  """Collect all the bits and pieces of an aligner analysis and compact them together into a
  single page, static HTML file with all the figures embedded as b64 encoded data"""
  open(htmlout, 'w').write(
    create_summary_page(tool, sample, graph, mismatcsv, summaryjson, novelplot, knownplot, mqplot, circleplot, matrixplot,
                        db_summary, known_d_plot, novel_d_plot)
  )


def create_summary_page(tool, sample, graph, mismatcsv, summaryjson, novelplot, knownplot, mqplot, circleplot, matrixplot,
                        db_summary, known_d_plot, novel_d_plot):
  """."""
  title = "{} - Aligner report".format(tool)

  html = ['<h3>Experiment details</h3>']
  html += metadata_table(tool, sample, graph)
  html += ['<h3>Accuracy summary</h3>']
  html += aligner_accuracy_summary_table(summaryjson)
  html += ['<h3>Plots</h3>']
  html += plots_table(novelplot, knownplot, known_d_plot, novel_d_plot, mqplot, circleplot, matrixplot)
  html += ['<h3>Alignment summary</h3>']
  html += misalignment_summary_table(mismatcsv)

  if db_summary is not None:
    html += ['<h3>Variant details</h3>']
    html += ['<pre>']
    html += [open(db_summary, 'r').read()]
    html += ['</pre>']

  return __html_template__.format(title, '\n'.join(html))


def metadata_table(tool, sample, graph):
  html = ['<table class="metadata">']
  html += ['<tr><td width=100px>Tool</td><td width=300px>{}</td></tr>'.format(tool)]
  html += ['<tr><td>Sample</td><td>{}</td></tr>'.format(sample)]
  html += ['<tr><td>Graph</td><td>{}</td></tr>'.format(graph)]
  html += ['</table>']
  return html


def aligner_accuracy_summary_table(summaryjson):
  """Load the data in the summary.json file and convert it into an html table
  "read_counts": {
    "known": {
      "DEL < 20 bp": 0,
      "INS >= 20 bp": 0,
      "DEL >= 20 bp": 0,
      "SNP": 0,
      "INS < 20 bp": 0
    },
    "novel": {
      "DEL < 20 bp": 25116,
      "INS >= 20 bp": 2714,
      "DEL >= 20 bp": 0,
      "SNP": 94043,
      "INS < 20 bp": 24105
    },
    "ref": 1120220
  },

  {
    "known": {
      "DEL < 20 bp": 0.0,
      "INS >= 20 bp": 0.0,
      "DEL >= 20 bp": 0.0,
      "SNP": 0.0,
      "INS < 20 bp": 0.0
    },
    "novel": {
      "DEL < 20 bp": 57.404610892982987,
      "INS >= 20 bp": 49.278152069297406,
      "DEL >= 20 bp": 0.0,
      "SNP": 89.022314423255935,
      "INS < 20 bp": 57.611608271398964
    },
    "ref": 93.741105267014063


  }"""
  data = json.load(open(summaryjson, 'r'))
  cols = ['DEL >= 20 bp', 'DEL < 20 bp', 'SNP', 'INS < 20 bp', 'INS >= 20 bp']

  html = []

  for k in sorted(data.keys()):
    if k == 'read_counts': continue

    this_data = data[k]

    html += ['<b>{}</b>'.format(k)]
    html += ['<table class="accuracy">']
    html += ['<tr><th width=100px></th><th width=60px>REF</th>'] + ['<th width=100px>{}</th>'.format(h) for h in cols] + ['</tr>']
    html += ['<tr><td>Novel</td><td rowspan=2>{:3.2f}</td>'.format(this_data['ref'])] + ['<td>{:3.2f}</td>'.format(this_data['novel'][c]) for c in cols] + ['</tr>']
    html += ['<tr><td>Known</td>'] + ['<td>{:3.2f}</td>'.format(this_data['known'][c]) for c in cols] + ['</tr>']
    html += ["</table>"]

  this_data = data['read_counts']

  html += ['<b>Read Counts</b>']
  html += ['<table class="accuracy">']
  html += ['<tr><th width=100px></th><th width=60px>REF</th>'] + ['<th width=100px>{}</th>'.format(h) for h in cols] + ['</tr>']
  html += ['<tr><td>Novel</td><td rowspan=2>{:,}</td>'.format(this_data['ref'])] + ['<td>{:,}</td>'.format(this_data['novel'][c]) for c in cols] + ['</tr>']
  html += ['<tr><td>Known</td>'] + ['<td>{:,}</td>'.format(this_data['known'][c]) for c in cols] + ['</tr>']
  html += ["</table>"]

  return html


def plots_table(novelplot, knownplot, known_d_plot, novel_d_plot, mqplot, circleplot, matrixplot):
  html = ["<table>"]
  html += ['<tr><th height=50px valign="bottom">{}</th></tr><tr><td style="text-align: center; vertical-align: middle;">{}</td></tr>'.format(fn[0], embed_image(fn[1]))
           for fn in [('Novel variants', novelplot), ('Known variants', knownplot),
                      ('Novel variants |d| plot', novel_d_plot), ('Known variants |d| plot', known_d_plot),
                      ('MQ', mqplot),
                      ('Misalignments plot (Circle)', circleplot), ('Misalignments plot (Matrix)', matrixplot)]
           if fn[1] is not None
           ]
  html += ["</table>"]
  return html


def misalignment_summary_table(mismatcsv):
  html = ['<table class="accuracy">']
  with open(mismatcsv, 'r') as fp:
    header = fp.readline().split(',')
    html += ['<tr>'] + ['<th style="transform:rotate(270deg); height:80px;width:50px;" rowspan=2>{}</th>'.format(h) for h in header[:2]]  # Source chrom and unmapped
    html += ['<th colspan={}>Mapped Chrom</th>'.format(len(header) - 2)]
    html += ['<th rowspan=2>Read Count</th>'] + ['</tr>']
    html += ['<tr>'] + ['<th>{}</th>'.format(h) for h in header[2:]] + ['</tr>']
    for line in fp:
      data = line.split(',')
      chrom = data[0]
      read_counts = [int(cnt) for cnt in data[1:]]
      total_reads = float(sum(read_counts))
      html += ['<tr>']
      html += ['<th>{}</th>'.format(chrom)] + ['<td>{:.3g} (%)</td>'.format(100 * cnt/total_reads) for cnt in read_counts]
      html += ['<td>{}</td>'.format(int(total_reads))]
      html += ['</tr>']
  html += ["</table>"]
  return html


def embed_image(fn):
  """Load image file and embed as base64."""
  if not fn: return ''
  with open(fn, 'rb') as fig_fp:
   encoded_string = base64.b64encode(fig_fp.read())
  return '<img src="data:image/png;base64,{}"></img>'.format(encoded_string)


if __name__ == '__main__':
  cli()