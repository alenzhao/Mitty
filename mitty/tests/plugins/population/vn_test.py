"""Run genomes generate and then analyze the resultant genomes file to make sure it satisfies the conditions we
set it"""
import tempfile
import os
import json

from click.testing import CliRunner

import mitty.lib.variants as vr
import mitty.genomes as genomes
import mitty.tests


def vn_test():
  """Test 'vn' population model"""
  _, param_file = tempfile.mkstemp(suffix='.json')
  _, db_file = tempfile.mkstemp(suffix='.hdf5')
  test_params = {
    "files": {
      "reference_file": mitty.tests.test_fasta_genome_file,
      "dbfile": db_file
    },
    "rng": {
      "master_seed": 12345
    },
    "population_model": {
      "vn": {
        "p_S": 0.104,
        "p_G": [0.1, 0.5, 0.9]
      }
    },
    "chromosomes": [1],
    "variant_models": [
      {
        "snp": {
          "p": 0.01
        }
      }
    ]
  }
  json.dump(test_params, open(param_file, 'w'))

  runner = CliRunner()
  result = runner.invoke(genomes.cli, ['generate', param_file])
  assert result.exit_code == 0, result
  assert os.path.exists(db_file)

  pop = vr.Population(fname=db_file, mode='r', in_memory=False)
  # ml = pop.get_variant_master_list(chrom=4)

  S = set(pop.get_sample_variant_index_for_chromosome(1, 'S')['index'])
  G0 = set(pop.get_sample_variant_index_for_chromosome(1, 'G0')['index'])
  G1 = set(pop.get_sample_variant_index_for_chromosome(1, 'G1')['index'])
  G2 = set(pop.get_sample_variant_index_for_chromosome(1, 'G2')['index'])

  known = S & G0

  assert known == S & G1
  assert known == S & G2

  novel = S - G0

  assert novel == S - G1
  assert novel == S - G2

  os.remove(param_file)
  os.remove(db_file)

# Need cleanup code to work even if test fails ...
# nosetests mitty.tests.plugins.population