from cpython cimport array
import array

cimport cython
import numpy as np
cimport numpy as np

cdef:
  int NOVEL_VARIANT = 1
  int KNOWN_VARIANT = 2


cdef advance_window(int[:] window,
                    np.int32_t[:] index,
                    int max_w_idx_1, int start, int stop,
                    np.int32_t[:] footprint_start,
                    np.int32_t[:] footprint_stop):
  while footprint_stop[index[window[0]]] < start:
    if window[0] < max_w_idx_1:
      window[0] += 1
    else:
      break

  while footprint_start[index[window[1]]] < stop:
    if window[1] < max_w_idx_1:
      window[1] += 1
    else:
      break
  # Since window is a list, it is changed in place


cdef find_variants_over_read(int[:] s_window, int[:] g_window,
                             int start, int stop, int r_chrom_cpy,
                             str[:] ref, str[:] alt,
                             np.int32_t[:] s_index, np.int8_t[:] s_gt,
                             np.int32_t[:] g_index,
                             np.int32_t[:] footprint_start, np.int32_t[:] footprint_stop,
                             bint has_graph):
  cdef int i
  return [(len(alt[s_index[i]]) - len(ref[s_index[i]]),
           NOVEL_VARIANT if not has_graph or s_index[i] not in g_index[g_window[0]:g_window[1] + 1] else KNOWN_VARIANT)
          for i in range(s_window[0], s_window[1])
          if (s_gt[i] == r_chrom_cpy or s_gt[i] == 2) and start < footprint_stop[s_index[i]] and stop >= footprint_start[s_index[i]]]


# Pure python
# 1014452 reads in 7.61s just to read/write
# 1014452 reads in 12.16s to also advance windows
# 1014452 reads in 15.59s to also compute read assignments
# 1014452 reads in 17.11s to do everything
#
# Cython + types
# 1014452 reads in 9.80s to do everythng
def read_assigner_iterator(in_bam, pop, chrom, sample_name, graph_name=None):
  """Setup features and prepare to iterate over the bam region. We expect in_bam to be an interator
  over reads from one chrom

  :param in_bam: a BAM region for given chromosome.
  :param pop: Population object
  :param chrom: [0, 1, 2, ...],
  :param sample_name: Name of sample
  :param graph_name: Name of graph. Leave None to indicate no graph
  :return:
  """
  cdef:
    bint has_graph = graph_name is not None
    int r_start, r_end, r_chrom_cpy, sample_size_1, graph_size_1
    int[:] sample_variant_window = array.array('i', [0, 0])
    int[:] graph_variant_window = array.array('i', [0, 0])
    np.int32_t[:] s_index, g_index, footprint_start, footprint_stop
    np.int8_t[:] s_gt
    str[:] ref, alt

  master_list = pop.get_variant_master_list(chrom=chrom).variants
  footprint_start, footprint_stop = master_list['pos'], master_list['stop']
  ref, alt = master_list['ref'], master_list['alt']

  sample = pop.get_sample_variant_index_for_chromosome(chrom=chrom, sample_name=sample_name)
  s_index, s_gt = sample['index'], sample['gt']
  sample_size_1 = sample.size - 1

  if graph_name is not None:
    graph = pop.get_sample_variant_index_for_chromosome(chrom=chrom, sample_name=graph_name)
    g_index = graph['index']
    graph_size_1 = graph.size - 1
  else:
    g_index, graph_size_1 = np.empty((0,), dtype='i4'), -1

  for r in in_bam:
    r_start = r.pos
    r_stop = r.get_tag('ZE')
    # rm_start = r.pnext
    # rm_stop = r.get_tag('Ze')
    r_chrom_cpy = r.get_tag('Zc')

    advance_window(sample_variant_window, s_index, sample_size_1, r_start, r_stop,
                   footprint_start, footprint_stop)
    if has_graph:
      advance_window(graph_variant_window, g_index, graph_size_1, r_start, r_stop,
                     footprint_start, footprint_stop)

    sample_vars = find_variants_over_read(
      sample_variant_window, graph_variant_window,
      r_start, r_stop, r_chrom_cpy,
      ref, alt,
      s_index, s_gt,
      g_index,
      footprint_start, footprint_stop,
      has_graph)

    if sample_vars:
      r.set_tag('Z0', [_v[0] for _v in sample_vars])
      r.set_tag('Z1', [_v[1] for _v in sample_vars])

    yield r