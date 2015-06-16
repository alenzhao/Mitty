#!python
from mitty.version import __version__

__cmd__ = """genomes ({:s}): Generate simulated genomes from a simulation parameter file.

Commandline::

  Usage:
    genomes generate <pfile>  [-v|-V] [-p]
    genomes dryrun <pfile>
    genomes inspect <dbfile>
    genomes inspect sfs <chrom> <dbfile>
    genomes write ((vcf|vcfs)|master) <dbfile> <out_prefix> [<serial>]... [-v|-V]
    genomes explain (parameters|(variant-model|spectrum-model|population-model) (all|<model_name>))
    genomes list (variant-model|spectrum-model|population-model)

  Options:
    generate                Create a database of genomes by running the models specified
    <pfile>                 Name for parameter file
    dryrun                  Don't run simulation, but print out useful info about it
    inspect                 Give summary information about a genome database
    <dbfile>                Name of genome database file
    sfs                     Print the site frequency spectrum
    <chrom>                 Which chromosome's site frequncy spectrum to print
    write                   Write out genome data from the database file in vcf format
    vcf                     Write out all genomes in one multi-sample vcf file
    vcfs                    Write out the genomes in separate single sample vcf files
    master                  Write  a VCF file with all the variants, no sample tag
    <serial>                Serial number of sample. Not needed and Ignored if 'master' is used. Leave blank with vcf/vcfs to write all samples
    <out_prefix>            Written files will have this prefix
    explain                 Explain the parameters/variant model/population model
    all                     Iterate over and explain all models
    <model_name>            Explain specific model
    list                    List the models
    -v                      Dump log messages
    -V                      Dump detailed log messages
    -p                      Show progress bar""".format(__version__)

__param__ = """Parameter file example::

  {
    # Path notes: an absolute path is left as is. A relative path is taken relative to the parameter file location
    "files": {
      "reference_dir": "/Users/kghose/Data/hg38/",  # If reference is chr1.fa, chr2.fa ... in this directory
      "reference_file": "/Users/kghose/Data/hg38/hg38.fa.gz",  # If reference is a single gz fasta file
      "dbfile": "Out/test.db"  # Output database file
    },
    "rng": {
      "master_seed": 1
    },
    "population_model": {
      "standard": {            # Name of sample chooser (population) model.
        "sample_size": 1,      # population model parameters
      }
    },
    "site_model": {
        "double_exp": {    # Name of model that handles the site frequency spectrum
          "k1": 0.1,       # spectrum model parameters
          "k2": 2.0,
          "p0": 0.001,
          "p1": 0.2,
          "bin_cnt": 30
        }
      }
    "chromosomes": [1, 2]  # Chromosomes to apply the models to
    "variant_models": [    # The list of variant models should come under this key
      {
        "snp": {           # name of the model.
          "p": 0.01        # Parameters required by the model
        }
      },
      {                    # We can chain as many models as we wish. We can also repeat models.
        "delete" : {
          "p": 0.01
        }
      }
    ]
  }"""
__doc__ = __cmd__ + __param__
# We split this up so that when we print the help, we can print it as separate pages. It got cluttered when printed all
# together. We want this to be the docstring because that is a nice format for our generated documentation

import json
import os
import time

import docopt
import numpy as np

import mitty.lib
import mitty.lib.util as mutil
import mitty.lib.io as mio
import mitty.lib.db as mdb
import mitty.lib.variants as vr

import logging
logger = logging.getLogger(__name__)


def cli():  # pragma: no cover
  """Serves as entry point for scripts"""
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__cmd__, ['-h'])
  else:
    cmd_args = docopt.docopt(__cmd__)  # Version string?

  level = logging.DEBUG if cmd_args['-V'] else logging.WARNING
  logging.basicConfig(level=level)

  if cmd_args['-v']:
    logger.setLevel(logging.DEBUG)

  if cmd_args['generate']:
    generate(cmd_args)
  elif cmd_args['dryrun']:
    dry_run(cmd_args)
  elif cmd_args['write']:
    write(cmd_args)
  elif cmd_args['explain']:
    explain(cmd_args)
  elif cmd_args['list']:
    print_list(cmd_args)
  elif cmd_args['sfs']:  # This should come before 'inspect', or we need to do an additional test
    print_sfs(cmd_args)
  elif cmd_args['inspect']:
    db_summary(cmd_args)


def dry_run(cmd_args):
  """Don't run simulation, but print out useful info about it"""
  params = json.load(open(cmd_args['<pfile>'], 'r'))
  sfs_model = load_site_frequency_model(params.get('site_model', None))
  if sfs_model:
    print('Site frequency model')
    print(sfs_model)
    p, f = sfs_model.get_spectrum()
    if len(p):  # Don't print this for a dud spectrum
      print('1/sum(p_i * f_i) = {:2.1f}'.format(1./(p*f).sum()))
    x, y = mutil.growth_curve_from_sfs(p, f)
    print('\nGrowth curve:\n')
    print('N\tv')
    for _x, _y in zip(x, y):
      print('{:d}\t{:f}'.format(_x, _y))
  else:
    print('No site model')


def generate(cmd_args):
  """Generate genomes based on the simulation parameter file

  :param cmd_args: from doc opt parsing
  """
  base_dir = os.path.dirname(cmd_args['<pfile>'])     # Other files will be with respect to this
  params = json.load(open(cmd_args['<pfile>'], 'r'))

  pop_db_name = mitty.lib.rpath(base_dir, params['files']['dbfile'])
  if os.path.exists(pop_db_name):
    logger.warning('Removed old simulation file')
    os.remove(pop_db_name)
  if not os.path.exists(os.path.dirname(pop_db_name)):
    logger.warning('Creating output directory {:s}'.format(pop_db_name))
    os.makedirs(os.path.dirname(pop_db_name))

  ref = mio.Fasta(multi_fasta=mitty.lib.rpath(base_dir, params['files'].get('reference_file', None)),
                  multi_dir=mitty.lib.rpath(base_dir, params['files'].get('reference_dir', None)))  # TODO: Ability to switch off persistence flag
  master_seed = int(params['rng']['master_seed'])
  assert 0 < master_seed < mitty.lib.SEED_MAX

  #sample_size = int(params['sample_size'])
  chromosomes = params['chromosomes']

  t0 = time.time()
  unique_variant_count, total_variant_count = \
    run_simulations(pop_db_name, ref, sfs_model=load_site_frequency_model(params.get('site_model', None)),
                    variant_models=load_variant_models(ref, params['variant_models']),
                    population_model=load_population_model(params.get('population_model', None), params),
                    chromosomes=chromosomes, master_seed=master_seed,
                    progress_bar_func=mitty.lib.progress_bar if cmd_args['-p'] else None)
  t1 = time.time()
  logger.debug('Took {:f}s'.format(t1 - t0))
  logger.debug('{:d} unique variants, {:d} variants in samples'.format(unique_variant_count, total_variant_count))


def run_simulations(pop_db_name, ref, sfs_model, variant_models=[], population_model=None,
                    chromosomes=[], master_seed=2, progress_bar_func=None):
  """Save the generated genome(s) into the database.

  :param pop_db_name:    name of database to save to
  :param ref:            Fasta object reference genome
  :param sfs_model:      Site frequency model object. If None, no re-balancing of variant probabilities will occur
  :param population_model population model object
  :param variant_models: list of variant model objects
  :param chromosomes:    list of chromosomes to simulate
  :param master_seed:    seed for all random number generations
  :param progress_bar_func: if a proper progress bar function is passed, this will show a progress bar as we complete
  """
  seed_rng = np.random.RandomState(seed=master_seed)
  conn = mdb.connect(db_name=pop_db_name)
  for n, meta in enumerate(ref.get_seq_metadata()):
    mdb.save_chromosome_metadata(conn, n + 1, **meta)

  p, f = sfs_model.get_spectrum() if sfs_model is not None else (None, None)
  unique_variant_count, total_variant_count = 0, 0
  for ch in chromosomes:
    ml = vr.VariantList()
    for m in variant_models:
      ml.add(*m.get_variants(ref=ref[ch]['seq'], chrom=ch, p=p, f=f, seed=seed_rng.randint(mutil.SEED_MAX)))
    ml.sort()
    if sfs_model is not None: ml.balance_probabilities(*sfs_model.get_spectrum())
    mdb.save_master_list(conn, ch, ml)
    unique_variant_count += len(ml)
    for gen, n, this_sample, frac_done in population_model.samples(chrom_no=ch, ml=ml, rng_seed=seed_rng.randint(mutil.SEED_MAX)):
      #this_sample = ml.generate_chromosome(rng)
      mdb.save_sample(conn, gen, n, ch, this_sample)
      total_variant_count += len(this_sample)
      if progress_bar_func is not None:
        progress_bar_func('Chrom {:d} '.format(ch), frac_done, 80)
    if progress_bar_func is not None: print('')
    #mdb.save_chromosome_metadata(conn, ch, ref[ch]['id'], len(ref[ch]['seq']), ref[ch]['md5'])
  conn.close()
  return unique_variant_count, total_variant_count


def load_site_frequency_model(sfs_model_json):
  if sfs_model_json is None:
    return None
  k, v = sfs_model_json.items()[0]
  return mitty.lib.load_sfs_plugin(k).Model(**v)


def load_variant_models(ref, model_param_json):
  """Given a json snippet corresponding to models and their parameters, load them"""
  return [mitty.lib.load_variant_plugin(k).Model(ref=ref, **v)
          for model_json in model_param_json
          for k, v in model_json.iteritems()]  # There really is only one key (the model name) and the value is the
                                               # parameter list


def load_population_model(pop_model_json, params={}):
  """Given a json snippet corresponding to the population model load it. If None, return the standard model

  :param pop_model_json: json snippet corresponding to the population model
  :param params:  the entire parameter json. For backward compatibility. If we have no pop model defined, as was normal
                  for Mitty v < 1.2.0, we default to the standard model (built-in). The single parameter, sample_size,
                  was defined in the main body of the parameter json, so we need that here.
  """
  if pop_model_json is None:
    pop_model_json = {'standard': {'sample_size': params['sample_size']}}  # params is assumed to have the key 'sample_size'
  k, v = pop_model_json.items()[0]
  return mitty.lib.load_pop_model_plugin(k).Model(**v)


def write(cmd_args):
  pop_db_name = cmd_args['<dbfile>']
  conn = mdb.connect(db_name=pop_db_name)
  out_prefix = cmd_args['<out_prefix>']
  if cmd_args['master']:
    write_master_vcf(conn, out_prefix)
  else:
    samples = cmd_args['<serial>']
    if len(samples) == 0:
      samples = [n for n in range(mdb.samples_in_db(conn))]
    if cmd_args['vcf']:
      write_samples_vcf(conn, samples, out_prefix)
    else:
      write_sample_vcfs(conn, samples, out_prefix)


def write_master_vcf(conn, out_prefix):
  contig_info = [mdb.load_chromosome_metadata(conn, ch[0])[1:] for ch in mdb.chromosomes_in_db(conn)]
  fname = out_prefix + '_master.vcf'
  with mio.vcf_for_writing(fname, [], contig_info) as fp:
    for ch in mdb.chromosomes_in_db(conn):
      ml = mdb.load_master_list(conn, ch[0])
      mio.write_chromosomes_to_vcf(fp, seq_id=ch[1], chrom_list=[], master_list=ml)


def write_sample_vcfs(conn, samples, out_prefix):
  contig_info = [mdb.load_chromosome_metadata(conn, ch[0])[1:] for ch in mdb.chromosomes_in_db(conn)]
  for spl in [int(s) for s in samples]:
    fname = out_prefix + '_s{:d}.vcf'.format(spl)
    with mio.vcf_for_writing(fname, ['s{:d}'.format(spl)], contig_info) as fp:
      for ch in mdb.chromosomes_in_db(conn):
        ml = mdb.load_master_list(conn, ch[0])
        mio.write_chromosomes_to_vcf(fp, seq_id=ch[1], chrom_list=[mdb.load_sample(conn, 0, spl, ch[0])], master_list=ml)


def write_samples_vcf(conn, samples, out_prefix):
  _samples = [int(s) for s in samples]
  contig_info = [mdb.load_chromosome_metadata(conn, ch[0])[1:] for ch in mdb.chromosomes_in_db(conn)]
  sample_header = ['s{:d}'.format(spl) for spl in _samples]
  fname = out_prefix + '_s{:d}.vcf'.format(spl)
  with mio.vcf_for_writing(fname, sample_header, contig_info) as fp:
    for ch in mdb.chromosomes_in_db(conn):
      ml = mdb.load_master_list(conn, ch[0])
      mio.write_chromosomes_to_vcf(fp, seq_id=ch[1], chrom_list=[mdb.load_sample(conn, 0, spl, ch[0]) for spl in _samples], master_list=ml)


def explain(cmd_args):
  if cmd_args['parameters']:
    print(__param__)
  else:
    kind = ['variant-model', 'spectrum-model', 'population-model']
    idx = [cmd_args[k] for k in kind].index(True)
    if cmd_args['all']:
      explain_all_models(kind=kind[idx])
    else:
      explain_model(name=cmd_args['<model_name>'], kind=kind[idx])


def explain_all_models(kind):
  discoverer = {
    'variant-model': mitty.lib.discover_all_variant_plugins,
    'spectrum-model': mitty.lib.discover_all_sfs_plugins,
    'population-model': mitty.lib.discover_all_pop_plugins
  }
  for name, mod_name in discoverer[kind]():
    explain_model(name, kind)


def explain_model(name, kind):
  """Load the given model of given kind and print description, parameter example and defaults

  :param name: name of model
  :param kind: 'variant', 'sfs' or 'pop'
  """
  loader = {
    'variant-model': mitty.lib.load_variant_plugin,
    'spectrum-model': mitty.lib.load_sfs_plugin,
    'population-model': mitty.lib.load_pop_model_plugin
  }
  try:
    mod = loader[kind](name)
  except ImportError as e:
    print('{0}: {1}'.format(name, e))
    print('Problem with loading {} model'.format(kind))
    return
  try:
    print('\n---- ' + name + ' (' + mod.__name__ + ') ----')
    print(mod._description)
    print(mitty.lib.model_init_signature_string(mod.Model.__init__))
  except AttributeError:
    print('No help for model "{:s}" available'.format(name))


def print_list(cmd_args):
  kind = ['variant-model', 'spectrum-model', 'population-model']
  idx = [cmd_args[k] for k in kind].index(True)
  print_model_list(kind[idx])


def print_model_list(kind):
  discoverer = {
    'variant-model': mitty.lib.discover_all_variant_plugins,
    'spectrum-model': mitty.lib.discover_all_sfs_plugins,
    'population-model': mitty.lib.discover_all_pop_plugins
  }
  print('\nAvailable {} models\n----------------'.format(kind))
  for name, mod_name in discoverer[kind]():
    print('- {:s} ({:s})\n'.format(name, mod_name))


def db_summary(cmd_args):
  """Print some useful information about the database

  :param cmd_args: parsed arguments
  """
  pop_db_name = cmd_args['<dbfile>']
  conn = mdb.connect(db_name=pop_db_name)
  chrom_list = mdb.chromosomes_in_db(conn)
  populated_chrom_list = mdb.populated_chromosomes_in_db(conn)
  n_s = mdb.samples_in_db(conn)
  variant_stats = np.empty((len(populated_chrom_list), 3), dtype=float)
  sample_max = 100
  for i, c in enumerate(populated_chrom_list):
    if n_s < sample_max:  # Take every sample
      ss = range(n_s)
    else:
      ss = np.random.randint(0, n_s, sample_max)
    s_len = np.empty(len(ss), dtype=float)
    for j, s in enumerate(ss):
      s_len[j] = len(mdb.load_sample(conn, 0, s, c[0]))
    _, _, seq_len, _ = mdb.load_chromosome_metadata(conn, c[0])
    variant_stats[i, :] = (s_len.mean(), s_len.std(), 1e6 * s_len.mean() / float(seq_len))

  print('{:s}'.format(pop_db_name))
  print('\tVariants in {:d} chromosomes (Genome has {:d})'.format(len(populated_chrom_list), len(chrom_list)))
  print('\t{:d} samples'.format(n_s))
  print('Unique variants in population')
  print('\tChrom\tVariants')
  for c in populated_chrom_list:
    print('\t{:d}\t{:d}'.format(c[0], mdb.variants_in_master_list(conn, c[0])))
  print('Variants in samples')
  print('\tChrom\tAvg variants\tStd variants\tVariants/megabase')
  for i, c in enumerate(populated_chrom_list):
    print('\t{:d}\t{:<9.2f}\t{:<9.2f}\t{:.1f}'.format(c[0], variant_stats[i, 0], variant_stats[i, 1], variant_stats[i, 2]))


def print_sfs(cmd_args):
  pop_db_name = cmd_args['<dbfile>']
  conn = mdb.connect(db_name=pop_db_name)
  chrom = int(cmd_args['<chrom>'])
  ml = mdb.load_master_list(conn, chrom)
  print('Site frequency spectrum for chrom {:d}'.format(chrom))
  print(ml)


if __name__ == "__main__":
  cli()