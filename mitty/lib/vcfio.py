import gzip
import time

import pysam
import numpy as np
import pandas as pd

import logging
logger = logging.getLogger(__name__)


_gt_dict = {
  '0|0': 0,
  '0|1': 1,
  '1|0': 2,
  '1|1': 3,
  '0/1': 4,  # We will alternate copy for this
  '1/1': 2
}

def _convert_gt(gt_string):
  """Function that converts the gt string (0|1, 1|0, 1|1) to a
  two bit genotype code. If the data is het un-phased it is
  alternately assigned a chromosome.

  :param row:
  :return:
  """
  gt = _gt_dict.get(gt_string[:3])  # We assume the GT value is the first entry in FORMAT
  if gt == 4:
    gt = _convert_gt.copy + 1
    _convert_gt.copy = not _convert_gt.copy
  return gt


def read_sample_from_vcf(fname, sample):
  gz = fname.endswith('gz')
  my_open = gzip.open if gz else open
  with my_open(fname, 'r') as f:
    for n, line in enumerate(f):
      if gz:
        line = line.decode()  # gzip returns bytes
      if line.startswith('#CHROM'):
        if sample not in line:
          raise RuntimeError('Sample {} does not exist in {}'.format(sample, fname))
        break
    else:
      raise RuntimeError('Malformed VCF {}'.format(fname))

  # TODO: handle large CSV files in chunks?
  return pd.read_csv(
    fname, delimiter='\t', header=n,
    usecols=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample],
    converters={sample: _convert_gt},
    dtype={'POS': np.uint32}).rename(columns={'#CHROM': 'CHROM', sample: 'GT'})


# This loads a whole file into memory.
# Right now, unphased variants always go into chrom copy 0|1
# TODO: Consider random assignment - but this really is a function of the truth VCF
def parse_variant_file(fname, sample):
  """Use pysam to read in a VCF file and create a recarray data structure suitable for read generation
  The salient features are
   - Separate arrays are created for each chromosome copy
   - The array contains not only the start reference base (traditional VCF POS) but also the
     last affected reference base. For SNP and INS this is the same as POS but for DEL this is the
     last base of the deletion. We need this when expanding sections of the genome

  :param fname:
  :param sample:
  :return: dict of numpy recarrays
  """
  mode = 'rb' if fname.endswith('bcf') else 'r'
  vcf_fp = pysam.VariantFile(fname, mode)
  vcf_fp.subset_samples([sample])
  return {
    k: parse_contig(vcf_fp, k)
    for k in vcf_fp.header.contigs.keys()
  }


# Hardcoded to diploid genomes for now.
def parse_contig(vcf_fp, contig_name):
  """

  :param vcf_fp:  pysam.VariantFile with .subset_samples applied, so we only get one sample
  :param contig_name:
  :return:
  """
  var_list = [v for v in vcf_fp.fetch(contig=contig_name)]
  return {
    gt: np.array(list(filter(None, (parse_variant(v, cpy) for v in var_list))),
             dtype=[('pos', 'u4'), ('stop', 'u4'), ('type', 'c'), ('len', 'u4'), ('alt', object)])
    for cpy, gt in zip([0, 1], ['1|0', '0|1'])
  }
  # The only way I could remember which was which was to explicitly write down 1|0 or 0|1


def parse_variant(v, cpy):
  """

  :param v:
  :return:
  """
  var = v.samples.values()[0]
  if var['GT'][cpy] == 0:  # Does not exist on this copy
    return None  # To be filtered out

  alt = var.alleles[cpy]
  if v.rlen > 1 and len(alt) > 1:
    raise ValueError('Complex variants can not be handled due to ambiguity in creating CIGARs. '
                     'Please filter out: {}:{} ref:{}, alt:{}'.format(v.contig, v.pos, v.ref, var.alleles))

  vtype, vlen = 'X', 0
  if v.rlen > 1:
    vtype, vlen, alt = 'D', v.rlen - 1, ''
  elif len(alt) > 1:
    vtype, vlen, alt = 'I', len(alt) - 1, alt[1:]

  return v.pos, v.stop, vtype, vlen, alt
  # v.stop is in 0 based indecies, but because it's exclusive we don't need to add 1


def prepare_variant_file(fname_in, sample, bed_fname, fname_out, write_mode='w'):
  """Prepare a variant file with only the given sample, complex variant calls filtered out, restricted to the given bed file

  :param fname_in:
  :param sample:
  :param bed_fname:
  :param fname_out:
  :return: - output is to file
  """
  def _complex_variant(_v):
    var = _v.samples.values()[0]
    for alt in var.alleles:
      if _v.rlen > 1 and len(alt) > 1:
        if _v.ref != alt:
          logger.debug('Filtered out {}:{} {} -> {}'.format(_v.contig, _v.pos, _v.ref, var.alleles))
          return True
    return False

  logger.debug('Starting filtering ...')
  t0 = time.time()

  mode = 'rb' if fname_in.endswith('bcf') else 'r'
  vcf_in = pysam.VariantFile(fname_in, mode)
  vcf_in.subset_samples([sample])
  vcf_out = pysam.VariantFile(fname_out, mode=write_mode, header=vcf_in.header)
  bed = list(map(lambda x: (x[0], int(x[1]), int(x[2])), map(lambda x: x.split(), open(bed_fname, 'r').readlines())))
  fltr_cnt = 0
  for region in bed:
    logger.debug('Filtering {}'.format(region))
    for v in vcf_in.fetch(contig=region[0], start=region[1], stop=region[2]):
      if _complex_variant(v):
        fltr_cnt += 1
        continue
      vcf_out.write(v)

  logger.debug('Filtered out {} complex variants'.format(fltr_cnt))
  t1 = time.time()
  logger.debug('Took {} s'.format(t1 - t0))