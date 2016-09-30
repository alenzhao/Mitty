"""Given a FASTQ file generate a table of Base Quality scores"""
from multiprocessing import Process, Queue
import logging

import numpy as np

import mitty.lib.fastq as fqi


logger = logging.getLogger(__name__)
__process_stop_code__ = 'SETECASTRONOMY'


# For debugging
def base_quality_single_threaded(fastq_fp, out_fp=None, f_size=None, max_reads=None):
  """Given a fastq file, read through the file

  :param fastq_fp:
  :param out_fp: If supplied, the score matrix will be written as a csv
  :param f_size: os.stat
  :param max_reads: Bug out after these many templates have been read
  :return:
  """
  #                 bp, phred
  score = np.zeros((500, 100), dtype=np.uint64)
  for template in fqi.read_fastq(fastq_fp, ipe=False, f_size=f_size, max_templates=max_reads):
    for n, bq in enumerate(template[0][3]):
      score[n, ord(bq) - 33] += 1

  if out_fp is not None:
   np.savetxt(out_fp, score, fmt='%d', delimiter=',')

  return score


def base_quality(fastq_fp, out_fp=None, threads=2, f_size=None, max_reads=None):
  """Given a fastq file, read through the file

  :param fastq_fp:
  :param out_fp: If supplied, the score matrix will be written as a csv
  :param threads: How many 'threads' to use
  :param f_size: os.stat
  :param max_reads: Bug out after these many templates have been read
  :return:
  """

  in_queue, out_queue = Queue(), Queue()

  # Start worker processes
  logger.debug('Starting {} threads'.format(threads))
  p_list = [Process(target=process_worker, args=(i, in_queue, out_queue)) for i in range(threads)]
  for p in p_list:
    p.start()

  # Burn through file
  logger.debug('Starting to read FASTQ file')
  for template in fqi.read_fastq(fastq_fp, ipe=False, f_size=f_size, max_templates=max_reads):
    in_queue.put(template)

  # Tell child processes to stop
  logger.debug('Stopping child processes')
  for i in range(threads):
    in_queue.put(__process_stop_code__)

  # Get results and add them
  score = out_queue.get()
  for i in range(threads - 1):
    score += out_queue.get()

  if out_fp is not None:
   np.savetxt(out_fp, score, fmt='%d', delimiter=',')

  # Wait for them to finish
  for p in p_list:
    p.join()

  return score


def process_worker(worker_no, in_queue, out_queue):
  """Create a bam fragment. This is designed to be a worker process for a multiprocessing
  pool, but can be tested without recourse to multiprocessing

  :param worker_no: an id for the work, not really used in computation
  :param in_queue:  an object with a get method that returns templates when called with next()
  :param out_queue: a queue to put the results on when done
  :return:
  """
  logger.debug('Worker {} starting'.format(worker_no))
  #                 bp, phred
  score = np.zeros((500, 100), dtype=np.uint64)
  for template in iter(in_queue.get, __process_stop_code__):
    for n, bq in enumerate(template[0][3]):
      score[n, ord(bq) - 33] += 1

  out_queue.put(score)
  logger.debug('Worker {} stopping'.format(worker_no))