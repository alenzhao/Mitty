import os


def rpath(base_dir, this_path):
  """Return this_path relative to base_dir, unless this_path is absolute"""
  return this_path if os.path.isabs(this_path) else os.path.normpath(os.path.join(base_dir, this_path))
