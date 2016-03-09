"""Stock hotspot plugin that distributes hotspots uniformly across a genome with a range of widths and heights.
The genome generator expects a list of tuples (pos, width, factor)

"""
import numpy as np


class Model:
  def __init__(
    self,
    n_hot=4,
    min_width=100,
    max_width=1000,
    min_height=2.0,
    max_height=10.0):
    self.n_hot = n_hot
    self.min_width = min_width
    self.max_width = max_width
    self.min_height = min_height
    self.max_height = max_height

  def hot_spots(self, ref, chrom, seed=1):
    rng = np.random.RandomState(seed=seed)
    ref_len = len(ref)
    return [
      (rng.uniform(0, ref_len),
       rng.uniform(low=self.min_width, high=self.max_width),
       rng.uniform(low=self.min_height, high=self.max_height))
      for _ in range(self.n_hot)
    ]