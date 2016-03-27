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
        text-align: left;
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
@click.argument('sample')
@click.argument('graph')
@click.argument('summaryjson', type=click.Path(exists=True))
@click.argument('novelplot', type=click.Path(exists=True))
@click.argument('knownplot', type=click.Path(exists=True))
@click.argument('mqplot', type=click.Path(exists=True))
@click.argument('htmlout', type=click.Path())
def cli(tool, sample, graph, summaryjson, novelplot, knownplot, mqplot, htmlout):
  """Collect all the bits and pieces of an aligner analysis and compact them together into a
  single page, static HTML file with all the figures embedded as b64 encoded data"""
  open(htmlout, 'w').write(
    create_summary_page(tool, sample, graph, summaryjson, novelplot, knownplot, mqplot)
  )


def create_summary_page(tool, sample, graph, summaryjson, novelplot, knownplot, mqplot):
  """."""
  title = "{} - Aligner report".format(tool)

  html = ['<h3>Experiment details</h3>']
  html += metadata_table(tool, sample, graph)
  html += ['<h3>Accuracy summary</h3>']
  html += aligner_accuracy_summary_table(summaryjson)
  html += ['<h3>Plots</h3>']
  html += plots_table(novelplot, knownplot, mqplot)

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

  html = ['<table class="accuracy">']
  html += ['<tr><th width=100px></th><th width=60px>REF</th>'] + ['<th width=100px>{}</th>'.format(h) for h in cols] + ['</tr>']
  html += ['<tr><td>Novel</td><td rowspan=2>{:3.2f}</td>'.format(data['ref'])] + ['<td>{:3.2f}</td>'.format(data['novel'][c]) for c in cols] + ['</tr>']
  html += ['<tr><td>Known</td>'] + ['<td>{:3.2f}</td>'.format(data['known'][c]) for c in cols] + ['</tr>']
  html += ["</table>"]

  return html


def plots_table(novelplot, knownplot, mqplot):
  html = ["<table>"]
  html += ['<tr><th height=50px valign="bottom">{}</th></tr><tr><td>{}</td></tr>'.format(fn[0], embed_image(fn[1]))
           for fn in [('Novel variants', novelplot), ('Known variants', knownplot), ('MQ', mqplot)]]
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