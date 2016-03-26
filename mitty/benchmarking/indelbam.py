"""Given a perfect BAM file from perfectbam use the information in the extended tags to work out
alignment accuracy parametrized by the sample indel lengths. Add indel information as extended tags
Return information in a json file"""
import time
import io
import cPickle
import json

import click
import pysam
import numpy as np

import mitty.version
import mitty.lib.variants as vr
import mitty.benchmarking.creed as creed

import logging
logger = logging.getLogger(__name__)


__extended_bam_tags_info__ = """
Z0  B    array of signed integer indicating variant(s) covered by this read.
         Value of integer indicates length of variant:
          0 = SNP,
          +N = N base insert
          -N = N base deletion
Z1  B    Array of variant flag(s)
         0x01 = 1 = variant found in sample only (novel variant)
         0x10 = 2 = variant found in graph and sample (known variant)
"""


def print_tags(ctx, param, value):
  if not value or ctx.resilient_parsing:
    return
  click.echo(__extended_bam_tags_info__)
  ctx.exit()


@click.command()
@click.argument('perbam', type=click.Path(exists=True))
@click.argument('gdb', type=click.Path(exists=True))
@click.argument('sample-name')
@click.argument('indelbam', type=click.Path())
@click.argument('indelpkl', type=click.Path())
@click.argument('summaryjson', type=click.Path())
@click.option('--indel-range', help='Maximum base pair count of indels we process', type=int, default=50)
@click.option('--graph-name', help='Name of graph used by aligner. Leave out for linear aligners')
@click.option('-p', is_flag=True, help='Show progress bar')
@click.option('-v', count=True, help='Verbosity level')
@click.option('--tags', is_flag=True, callback=print_tags, expose_value=False, is_eager=True, help='Print documentation for extended BAM tags')
def cli(perbam, gdb, indelbam, sample_name, indelpkl, summaryjson, indel_range, graph_name, p, v):
  """For each read in a BAM, use extended tags to indicate which indel categories it belongs to.
  Tags are added to original BAM in place"""
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  perbam_fp = pysam.AlignmentFile(perbam, 'rb')
  outbam_fp = pysam.AlignmentFile(indelbam, 'wb', header=create_bam_header(perbam_fp))

  total_read_count = perbam_fp.mapped + perbam_fp.unmapped  # Sadly, this is only approximate
  progress_bar_update_interval = int(0.01 * total_read_count)

  t0 = time.time()
  with click.progressbar(
    length=total_read_count, label='Analyzing indels', file=None if p else io.BytesIO()) as bar:
    pop = vr.Population(gdb)
    read_counter = creed.init_read_counts(indel_range)
    for ch in pop.get_chromosome_list():
      f_chrom_id = perbam_fp.header['SQ'][ch - 1]['SN']
      assign_features_to_reads(
        in_bam=perbam_fp.fetch(region=f_chrom_id), out_bam=outbam_fp,
        pop=pop, chrom=ch, read_counter=read_counter,
        sample_name=sample_name, graph_name=graph_name,
        progress_callback=bar.update,
        progress_bar_update_interval=progress_bar_update_interval)

    alignment_data = creed.read_counts_dereference(read_counter)

    cPickle.dump(alignment_data, open(indelpkl, 'wb'), protocol=cPickle.HIGHEST_PROTOCOL)
    json.dump(summarize(alignment_data), open(summaryjson, 'w'), indent=2)

  t1 = time.time()
  logger.debug('Analyzed {:d} reads in {:2.2f}s'.format(total_read_count, t1 - t0))

  t0 = time.time()
  pysam.index(indelbam)
  t1 = time.time()
  logger.debug('Indexed BAM in {:2.2f}s'.format(t1 - t0))


def create_bam_header(perbam_fp):
  new_header = perbam_fp.header
  pp = new_header['PG'][-1]['ID'] if 'PG' in new_header else None
  if 'PG' not in new_header:
    new_header['PG'] = []
  new_header['PG'].append({
    'CL': 'indelbam ....',
    'ID': 'mitty-indelbam',
    'PN': 'indelbam',
    'VN': mitty.version.__version__,
    'PP': pp,
    'DS': 'Attach extended tags describing variants covering reads'
  })
  return new_header


def assign_features_to_reads(
  in_bam, out_bam, pop, chrom, read_counter,
  sample_name, graph_name=None,
  progress_callback=None, progress_bar_update_interval=None):
  """

  :param in_bam: Iterator over input bam file
  :param out_bam: Output BAM file
  :param pop: Populaion object
  :param chrom: [0, 1, 2, ...],
  :param read_counter:
  :param sample_name: Name of sample
  :param graph_name: Name of graph. Leave None to indicate no graph
  :param progress_callback: if supplied, this function is called periodically with a single
                            argument, which is an int indicating how many reads have been
                            processed since the last call to this function
  :param progress_bar_update_interval: How many reads to process before updating prog bar
  :return:

  We assume all data here - the reads and the start arrays are sorted by position. This allows us
  to do a window based comparison that is more efficient than the n^2 computation un sorted data
  would entail
  """
  for cnt, r in enumerate(creed.read_assigner_iterator(in_bam, pop, chrom,
                                                       read_counter,
                                                       sample_name, graph_name)):
    out_bam.write(r)
    if (cnt + 1) % progress_bar_update_interval == 0:
      progress_callback(progress_bar_update_interval)
  progress_callback(progress_bar_update_interval)


indel_categories = [('DEL >= 20 bp', -100000, -20), ('DEL < 20 bp', -20, -1),
                    ('INS < 20 bp', 0, 19), ('INS >= 20 bp', 19, 100000)]


def summarize(alignment_data):
  """Simply tally up percentage corrects for Ref, SNP and Indel"""
  ad = alignment_data
  ref_pc = ad['ref_r_cnt'][0] / float(ad['ref_r_cnt'].sum()) * 100

  v_len = ad['v_len']
  idx_snp = np.where(v_len == 0)[0][0]

  known_v_cnt = ad['known_v_r_cnt']
  novel_v_cnt = ad['novel_v_r_cnt']

  known_snp_pc = known_v_cnt[idx_snp][0] / float(known_v_cnt[idx_snp].sum() or 1) * 100
  novel_snp_pc = novel_v_cnt[idx_snp][0] / float(novel_v_cnt[idx_snp].sum() or 1) * 100

  stats = {
    'ref': ref_pc,
    'known': {
      'SNP': known_snp_pc
    },
    'novel': {
      'SNP': novel_snp_pc
    }
  }
  for k, r0, r1 in indel_categories:
    indel_idx = np.where((r0 < v_len) & (v_len < r1))[0]
    stats['known'][k] = known_v_cnt[indel_idx, 0].sum() / float(known_v_cnt[indel_idx, :].sum() or 1) * 100
    stats['novel'][k] = novel_v_cnt[indel_idx, 0].sum() / float(novel_v_cnt[indel_idx, :].sum() or 1) * 100

  return stats


if __name__ == '__main__':
  cli()