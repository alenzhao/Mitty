"""This module contains functions that generate simulated reads. The module can be called as a script as well. This is
useful for creating test data for MGR algorithms/data formats

Commandline::

  Usage:
    reads [localreads]  --wg=WG  --out=OUT  --paramfile=PFILE  [--corrupt]  [--fastq] [--read_border_size=RBS] [--reads_per_block=BL] [--master_seed=MS] [-v]

  Options:
    localreads              Take reads around insertions only
    --wg=WG                 Whole genome file
    --out=OUT               Output file name prefix
    --read_border_size=RBS  How many bases before and after the insertion do we include in our window [default: 500]
    --paramfile=PFILE       Name for parameter file
    --corrupt               Write out corrupted reads too.
    --fastq                 Write as FASTQ instead of BAM (simulated_reads.fastq)
    --reads_per_block=BL    Generate these many reads at a time (Adjust to machine resources). [default: 100000]
    --master_seed=MS        If this is specified this passes a master seed to the read plugin.
                            This overrides any individual seeds specified by the parameter file.
    -v                      Dump detailed logger messages

Parameter file example::

  {
    "take reads from": [1,2],
    "coverage": 5.0,
    "output_file_prefix": "sim_reads",
    "read_model": "simple_reads",
    "model_params": {
      "paired": false,
      "read_len": 100,
      "template_len": 250
    }
  }



1. The quality scores are in Phred scale (as specified in the SAM spec)
2. We supply the prefix of output file name in the parameter file . Say we set this as sim_reads.
   The perfect reads will be saved to sim_reads.bam (or sim_reads.fastq). If we ask for corrupted reads
   we will get the corrupted reads in the file sim_reads_c.fastq.
   A text sidecar file sim_reads.info will always be saved with simulation parameters.

3 *** Can probably refactor the whole file for better readability **
  *** see if we can write shorter functions, see if we can reduce the number of parameters/bundle them ***
  *** further efficiency gains will probably be minimal - the bottle neck function has been cythonized ***

"""

__explain__ = """
Example parameter file .json

{
    "whole genome file": "chimera.h5",
    "take reads from": [1,2],
    "coverage": 5.0,
    "output_file_prefix": "sim_reads",
    "read_model": "simple_reads",
    "model_params": {
        "paired": false,
        "read_len": 100,
        "template_len": 250
    }
}

* If "take reads from" is set to none reads will be generated from all chromosomes. Otherwise reads will be taken only
  from the specified chromosomes
* If the pos data for a particular chromosome is not present, we will assume that the sequence is the reference sequence

Qname "cheat" string for reads

The qname of an unpaired read is written as

    ch:cp|rN|POS|CIGAR|PABS

while that of paired reads are written as

    ch:cp|rN|POS1|CIGAR1|PABS1|POS2|CIGAR2|PABS2

Where

    ch    - chromosome number
    cp    - chromosome copy
    N     - number of the read,
    POS   - correct position of read
    CIGAR - correct CIGAR
    PABS  - actual position of read from physical sequence. For reference sequences this will match POS,
            for variants, this will not
"""

__version__ = '0.4.0'

import numpy
import pyximport; pyximport.install(setup_args={"include_dirs": numpy.get_include()})
from roll_cigar import roll_cigar
import h5py
import string
import os
import imp
import json
import docopt
import pysam  # Needed to write BAM files
import logging
logger = logging.getLogger(__name__)

DNA_complement = string.maketrans('ATCGN', 'TAGCN')


#TODO: change this to work on list of reads so we avoid function overhead
def interpret_read_qname(read):
  """Given a read generated by this tool give us the information about the perfect alignment stored in
  that tool.
  The qname looks like ->  ch:cp|rN|POS1|CIGAR1|PABS1|POS2|CIGAR2|PABS2
  """
  cheat_answer = read.qname.split('|')
  correct_chrom_no, correct_chrom_copy = [int(k) for k in cheat_answer[0].split(':')]
  correct_pos = int(cheat_answer[2])
  correct_cigar = cheat_answer[3]
  absolute_pos = int(cheat_answer[4])
  if read.flag & 0x01:  # Paired reads
    if read.flag & 0x80:  # Second end (mate)
      correct_pos = int(cheat_answer[5])
      correct_cigar = cheat_answer[6]
      absolute_pos = int(cheat_answer[7])
  return correct_chrom_no, correct_chrom_copy, correct_pos, correct_cigar, absolute_pos


def open_reads_files(out_prefix, corrupted_reads=False, save_as_bam=True):
  """Open output files (fastq or bam) and write headers as needed."""
  file_handles = {}
  if save_as_bam:  # BAM
    perfect_reads_fname = out_prefix + '.bam'
    bam_hdr = {'HD': {'VN': '1.4'},
               'SQ': [{'LN': 1, 'SN': 'raw reads',
                       'SP': 'reads.py {:s}'.format(__version__)}]}
    file_handles['perfect'] = pysam.Samfile(perfect_reads_fname, 'wb', header=bam_hdr)  # Write binary BAM with header
    if corrupted_reads:
      corrupted_reads_fname = out_prefix + '_c.bam'
      bam_hdr['SQ'] = [{'LN': 1, 'SN': 'raw corrupted reads',
                        'SP': 'reads.py {:s}'.format(__version__)}]
      file_handles['corrupted'] = \
        pysam.Samfile(corrupted_reads_fname, 'wb', header=bam_hdr)  # Write binary BAM with header
  else:  # FASTQ
    perfect_reads_fname = out_prefix + '.fastq'
    file_handles['perfect'] = open(perfect_reads_fname, 'w')  # File handles for FASTQ
    if corrupted_reads:
      corrupted_reads_fname = out_prefix + '_c.fastq'
      file_handles['corrupted'] = open(corrupted_reads_fname, 'w')

  file_handles['bam_file?'] = save_as_bam
  logger.debug('Saving perfect reads to {:s}'.format(perfect_reads_fname))
  try:
    logger.debug('Saving corrupted reads to {:s}'.format(corrupted_reads_fname))
  except NameError:
    pass
  return file_handles


def close_reads_files(file_handles):
  file_handles['perfect'].close()
  try:
    file_handles['corrupted'].close()
  except KeyError:
    pass
  #for k, v in file_handles.iteritems():
  #  v.close()


# ## This function is the bottleneck, taking 1ms to run per call
# def roll_cigar(this_read, p_arr):
#   """
#   You can 'read along' with these tests from the Readme developers section
#
#   Test a fully matching read
#   >>> t_read = ['ATTG','~~~~', 0]; \
#   p_arr = [1, 2, 3, 4, 5, 5, 5, 6, 9]; \
#   roll_cigar(t_read, p_arr)
#   (1, '4M')
#
#   Test for read with insert
#   >>> t_read = ['TTGT', '~~~~', 1]; \
#   roll_cigar(t_read, p_arr)
#   (2, '3M1I')
#
#   Another test for read with insert
#   >>> t_read = ['TGTT', '~~~~', 2]; \
#   roll_cigar(t_read, p_arr)
#   (3, '2M2I')
#
#   Test for read with delete at end - should not show up in CIGAR
#   >>> t_read = ['TTAC', '~~~~', 4]; \
#   roll_cigar(t_read, p_arr)
#   (5, '2I2M')
#
#   Test for read spanning a deletion - should get a delete
#   >>> t_read = ['ACAC', '~~~~', 0]; \
#   p_arr = [1, 2, 5, 6, 7, 8, 9]; \
#   roll_cigar(t_read, p_arr)
#   (1, '2M2D2M')
#
#   We actually missed this case: read with one matching base and then a delete
#   >>> t_read = ['CACT', '~~~~', 1]; \
#   roll_cigar(t_read, p_arr)
#   (2, '1M2D3M')
#
#   Test for an unmapped read: pos and cigars should be None
#   >>> t_read = ['AATT', '~~~~', 2]; \
#   p_arr = [1, 2, 3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9]; \
#   roll_cigar(t_read, p_arr)
#   (0, '')
#
#   """
#   mapped = False
#   cigar = ''
#   coord = this_read[2]
#   counter = 0
#   cigar_fragment = None
#   for n in range(coord, coord + len(this_read[0])):
#     dp = p_arr[n+1] - p_arr[n]
#     if dp == 1:
#       mapped = True  # As long as we have one M we are a mapped read
#       if cigar_fragment != 'M':
#         if counter > 0:  # Flush
#           cigar += '{:d}{:s}'.format(counter, cigar_fragment)
#           counter = 0
#       cigar_fragment = 'M'
#       counter += 1
#     elif dp == 0:
#       if cigar_fragment != 'I':
#         if counter > 0:  # Flush
#           cigar += '{:d}{:s}'.format(counter, cigar_fragment)
#           counter = 0
#       cigar_fragment = 'I'
#       counter += 1
#     elif dp > 1:
#       mapped = True  # As long as we have one M we are a mapped read
#       if cigar_fragment != 'M':
#         if counter > 0:  # Flush
#           cigar += '{:d}{:s}'.format(counter, cigar_fragment)
#           counter = 0
#       cigar_fragment = 'M'  # We need to set this because we could be at the start of a read and type = None still
#       counter += 1
#       cigar += '{:d}{:s}'.format(counter, cigar_fragment)
#       cigar_fragment = 'D'
#       counter = dp - 1
#
#   if cigar_fragment != 'D':  # Flush all but 'D'. We only write D if we cross a D boundary
#     cigar += '{:d}{:s}'.format(counter, cigar_fragment)
#
#   if mapped:
#     align_pos = p_arr[coord]
#   else:
#     align_pos = 0
#     cigar = ''
#
#   return align_pos, cigar


# The qname of an unpaired read is written as
# ch:cp|rN|POS1|PABS1|CIGAR1|
# while that of paired reads are written as
# ch:cp|rN|POS1|PABS1|CIGAR1|POS2|CIGAR2|PABS2
# This function fills out the POS and CIGAR in the qname
def roll(these_reads, pos_array):
  """Given a list of reads return us the reads with the coordinate replaced by a 'POS:CIGAR' string.
  CIGAR, create the CIGAR, roll, get it? Oh, alright, I tried.

  Inputs:
                                 _________ ( seq_str, quality_str, coordinate)
    reads     -  [              /
                  [[ ... ], [ ... ]],
                  [[ ... ], [ ... ]], -> inner list = 2 elements if paired reads, 1 otherwise
                       .
                       .
                       .
                 ] -> outer list = number of reads

    pos_array -  1 d array

  Output:
    No output - changes reads in place
  """
  def roll_null(these_reads):
    """Simple convenience function - for if we have null reads ."""
    paired = True if len(these_reads[0]) == 2 else False
    for this_read in these_reads:
      qname = '{:d}|{:d}M|{:d}'.format(this_read[0][2] + 1, len(this_read[0][0]), this_read[0][2] + 1)
      if paired:
        qname += '|{:d}|{:d}M|{:d}'.format(this_read[1][2] + 1, len(this_read[1][0]), this_read[1][2] + 1)
        this_read[1][2] = qname
      this_read[0][2] = qname

  if len(these_reads) == 0: return
  paired = len(these_reads[0]) == 2
  # If these reads are from the reference sequence, then return CIGARs from the null read
  if pos_array is None:
    return roll_null(these_reads)

  for this_read in these_reads:
    pos, cigar = roll_cigar(this_read[0], pos_array)
    qname = '{:d}|{:s}|{:d}'.format(pos, cigar, this_read[0][2] + 1)
    if paired:
      pos, cigar = roll_cigar(this_read[1], pos_array)
      qname += '|{:d}|{:s}|{:d}'.format(pos, cigar, this_read[1][2] + 1)
      this_read[1][2] = qname
    this_read[0][2] = qname


def save_reads(fh, these_reads, chrom_key, offset, save_as_bam=True):
  if save_as_bam:
    save_reads_to_bam(fh, these_reads, chrom_key, offset)
  else:
    save_reads_to_fastq(fh, these_reads, chrom_key, offset)


# Didn't refactor these functions as that would involve another function call within a loop
#TODO: Implement option to write short qname
def save_reads_to_bam(bam_file, these_reads, chrom_key, first_read_serial):
  """
  Inputs:
    bam_file  -  bam file handle

                                 _________ ( seq_str, quality_str, 'POS:CIGAR ...' string to put in qname)
    reads     -  [              /
                  [( ... ), ( ...)],
                  [( ... ), ( ...)], -> inner list = 2 elements if paired reads, 1 otherwise
                       .
                       .
                       .
                 ] -> outer list = number of reads
    chrom_key  - chrom key string used to get sequence from .wg file e.g. '1:2'
    first_read_serial  - serial number of first read in this list. Needed since we give each template a unique name
  """
  if len(these_reads) == 0: return
  paired = True if len(these_reads[0]) == 2 else False
  for ser_no, this_read in enumerate(these_reads):
    ar = pysam.AlignedRead()
    ar.qname = '{:s}|r{:d}|{:s}'.format(chrom_key, first_read_serial + ser_no, this_read[0][2])
    ar.seq = this_read[0][0]
    ar.qual = this_read[0][1]
    if paired:
      ar.flag = 0x41  # end1 0x01 flag has to be set to indicate multiple segments
    bam_file.write(ar)

    if paired:
      ar.seq = this_read[1][0]
      ar.qual = this_read[1][1]
      ar.flag = 0x81  # end2 0x01 flag has to be set to indicate multiple segments
      bam_file.write(ar)


def save_reads_to_fastq(fastq_file_handle, these_reads, chrom_key, first_read_serial):
  """Given a list of sequences and their quality write the read data to an already opened text file. This saves data
  in interleaved form if paired reads are present

    bam_file  -  bam file handle

                                 _________ ( seq_str, quality_str, 'POS:CIGAR ...' string to put in qname)
    reads     -  [              /
                  [( ... ), ( ...)],
                  [( ... ), ( ...)], -> inner list = 2 elements if paired reads, 1 otherwise
                       .
                       .
                       .
                 ] -> outer list = number of reads
    chrom_key  - chrom key string used to get sequence from .wg file e.g. '1:2'
    first_read_serial  - serial number of first read in this list. Needed since we give each template a unique name
  """
  if len(these_reads) == 0: return
  paired = True if len(these_reads[0]) == 2 else False
  for ser_no, this_read in enumerate(these_reads):
    qname = '{:s}|r{:d}|{:s}'.format(chrom_key, first_read_serial + ser_no, this_read[0][2])
    seq = this_read[0][0]
    qual = this_read[0][1]
    fastq_file_handle.write('@{:s}\n{:s}\n+\n{:s}\n'.format(qname, seq, qual))
    if paired:
      seq = this_read[1][0]
      qual = this_read[1][1]
      fastq_file_handle.write('@{:s}\n{:s}\n+\n{:s}\n'.format(qname, seq, qual))


def add_reads_to_file(chrom=1,
                      cpy=1,
                      seq=[None,None],
                      seq_pos=None,
                      read_start=0,
                      read_stop=None,
                      num_reads=1000,
                      reads_per_call=100000,
                      generate_corrupt_reads=False,
                      read_model=None,
                      model_params={},
                      file_handles={}):

  reads = read_model.read_generator(seq=seq,
                                    read_start=read_start,
                                    read_stop=read_stop,
                                    num_reads=num_reads,
                                    reads_per_call=reads_per_call,
                                    generate_corrupt_reads=generate_corrupt_reads,
                                    **model_params)
  logger.debug('Generating reads')
  chrom_key = '{:d}:{:d}'.format(chrom, cpy)
  current_template_count = 0
  for perfect_reads, corrupt_reads in reads:
    roll(perfect_reads, seq_pos)  # The function modifies reads in place to fill out the POS and CIGAR
    logger.debug('Computed POS and CIGARs')
    save_reads(file_handles['perfect'], perfect_reads, chrom_key, current_template_count, file_handles['bam_file?'])
    #save_reads needs current_template_count because we put in a read serial number in the qname
    logger.debug('Saved perfect reads')
    if generate_corrupt_reads:
      for c_template, p_template in zip(corrupt_reads, perfect_reads):
        for c_read in c_template:
          c_read[2] = p_template[0][2]  # Copy over the qname
      save_reads(file_handles['corrupted'], corrupt_reads, chrom_key, current_template_count, file_handles['bam_file?'])
      logger.debug('Saved corrupted reads')
    current_template_count += len(perfect_reads)
    generated_reads = current_template_count * len(perfect_reads[0])  # This last gives us reads per template
    logger.debug('{:d} reads of {:d} done ({:d}%)'.
                 format(generated_reads, num_reads, int(100 * generated_reads / float(num_reads))))


def get_sequence_and_complement(ref_fp, chrom, cpy):
  logger.debug('Loading sequence')
  seq = ref_fp['sequence/{:d}/{:d}'.format(chrom, cpy)][:].tostring()
  logger.debug('Computing complement')
  complement_seq = seq.translate(DNA_complement)
  seq_len = len(seq)
  try:
    seq_pos = ref_fp['pos/{:d}/{:d}'.format(chrom, cpy)][:]
  except KeyError:
    logger.debug('Chrom {:d}:{:d} is same as reference'.format(chrom, cpy))
    seq_pos = None
  return seq, complement_seq, seq_len, seq_pos


def reads_around_insertions(args, params, read_model):
  """Generate reads only around insertions.

  We take the list of chromosomes
  chromosomes = [(chrom, (cpy,cpy ...)) ...] telling us which chromosomes and copies to do
  Numbering starts from 1

  """
  save_as_bam = not args['--fastq']
  write_corrupted = args['--corrupt']  # If True, corrupted reads will be written out
  read_border_size = int(args['--read_border_size'])
  coverage = params['coverage']
  output_file_prefix = args['--out']
  model_params = params['model_params']
  with h5py.File(args['--wg'], 'r') as ref_fp:
    chrom_list = params["take reads from"] or ref_fp['sequence'].keys()
    for chrom in chrom_list:
      if str(chrom) not in ref_fp['sequence'].keys():
        logger.warning('Chromosome #{:d} not in file'.format(chrom))
        continue
      logger.debug('Generating reads from chromosome {:d}'.format(chrom))
      try:
        vp = ref_fp['variants/pos/{:d}'.format(chrom)]
      except KeyError:
        vp = None

      for cpy in [int(c) for c in ref_fp['sequence/{:d}'.format(chrom)].keys()]:
        seq, complement_seq, seq_len, seq_pos = get_sequence_and_complement(ref_fp, chrom, cpy)
        if vp is None:
          logger.warning('Local reads around variants requested, but no variants detected')
          continue
        tvp = vp[:, cpy - 1, :]
        insertion_count = tvp.shape[0]
        for n in range(insertion_count):
          logger.debug('Generating reads from insertion {:d} of {:d}, len {:d}'.format(n, insertion_count, tvp[n, 1] - tvp[n, 0]))
          read_start = max(0, tvp[n, 0] - read_border_size)
          read_stop = min(seq_len, tvp[n, 1] + read_border_size)
          total_reads = int((read_stop - read_start) * coverage / read_model.average_read_len(**model_params))
          this_output_file_prefix = '{:s}_chrom{:d}_cpy{:d}_vn{:d}'.format(output_file_prefix, chrom, int(cpy), n)
          reads_file_handles = open_reads_files(this_output_file_prefix, write_corrupted, save_as_bam)
          add_reads_to_file(chrom=chrom, cpy=int(cpy), seq=[seq, complement_seq], seq_pos=seq_pos, read_start=read_start, read_stop=read_stop,
                            num_reads=total_reads,
                            reads_per_call=int(args['--reads_per_block']),
                            generate_corrupt_reads=write_corrupted,
                            read_model=read_model,
                            model_params=model_params,
                            file_handles=reads_file_handles)
          close_reads_files(reads_file_handles)


def whole_genome_reads(args, params, read_model):
  """Generate reads from the whole genome."""
  # Generate a dictionary of file handles for perfect and corrupted reads (if needed)
  save_as_bam = not args['--fastq']
  write_corrupted = args['--corrupt']  # If True, corrupted reads will be written out
  reads_file_handles = open_reads_files(args['--out'], write_corrupted, save_as_bam)
  coverage = params['coverage']
  model_params = params['model_params']
  with h5py.File(args['--wg'], 'r') as ref_fp:
    chrom_list = params["take reads from"] or ref_fp['sequence'].keys()
    for chrom in chrom_list:
      if str(chrom) not in ref_fp['sequence'].keys():
        logger.warning('Chromosome #{:d} not in file'.format(chrom))
        continue
      logger.debug('Generating reads from chromosome {:d}'.format(chrom))
      for cpy in ref_fp['sequence/{:d}'.format(chrom)].keys():
        seq, complement_seq, seq_len, seq_pos = get_sequence_and_complement(ref_fp, chrom, int(cpy))
        read_start = 0
        read_stop = seq_len
        total_reads = int((read_stop - read_start) * coverage / read_model.average_read_len(**model_params))
        add_reads_to_file(chrom=chrom, cpy=int(cpy), seq=[seq, complement_seq], seq_pos=seq_pos, read_start=read_start, read_stop=read_stop,
                          num_reads=total_reads,
                          reads_per_call=int(args['--reads_per_block']),
                          generate_corrupt_reads=write_corrupted,
                          read_model=read_model,
                          model_params=model_params,
                          file_handles=reads_file_handles)
  close_reads_files(reads_file_handles)


def main(args):
  # Load parameter file
  params = json.load(open(args['--paramfile'], 'r'))

  # Load the read model from the plugins directory
  plugin_dir = os.path.join(os.path.dirname(__file__), 'Plugins', 'Reads')
  fp, pathname, description = imp.find_module(params['read_model'] + '_plugin', [plugin_dir])
  read_model = imp.load_module('readmodel', fp, pathname, description)

  if args['--master_seed']:
    params['model_params']['master_seed'] = int(args['--master_seed'])

  if args['localreads']:
    reads_around_insertions(args, params, read_model)
  else:
    whole_genome_reads(args, params, read_model)

  with open(args['--out'] + '.info', 'w') as f:
    f.write('Command line\n-------------\n')
    f.write(json.dumps(args, indent=4))
    f.write('\n\nParameters\n------------\n')
    f.write(json.dumps(params, indent=4))


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    cmd_args = docopt.docopt(__doc__, version=__version__)

  level = logging.DEBUG if cmd_args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  main(cmd_args)
