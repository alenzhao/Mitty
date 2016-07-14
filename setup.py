from setuptools import setup, find_packages

__version__ = eval(open('mitty/version.py').read().split('=')[1])
setup(
    name='mitty',
    version=__version__,
    description='Simulator for genomic data',
    author='Seven Bridges Genomics',
    author_email='kaushik.ghose@sbgenomics.com',
    packages=find_packages(include=['mitty*']),
    include_package_data=True,
    entry_points={
      # Register the built in plugins
      'mitty.plugins.sfs': ['double_exp = mitty.plugins.site_frequency.double_exp'],
      'mitty.plugins.variants': ['snp = mitty.plugins.variants.snp_plugin',
                                 'delete = mitty.plugins.variants.delete_plugin',
                                 'uniformdel = mitty.plugins.variants.uniform_deletions',
                                 'uniformins = mitty.plugins.variants.uniform_insertions',
                                 'insert = mitty.plugins.variants.insert_plugin',
                                 #'inversion = mitty.plugins.variants.inversion_plugin',
                                 #'low_entropy_insert = mitty.plugins.variants.low_entropy_insert_plugin'
                                 ],
      'mitty.plugins.variants.hotspot': ['uniform = mitty.plugins.variants.hotspots.uniformhot'],
      'mitty.plugins.population': ['standard = mitty.plugins.population.standard',
                                   'vn = mitty.plugins.population.vn'],
      'mitty.plugins.reads': ['simple_sequential = mitty.plugins.reads.simple_sequential_plugin',
                              'simple_se = mitty.plugins.reads.simple_se',
                              'simple_illumina = mitty.plugins.reads.simple_illumina_plugin'],
      # Command line scripts
      'console_scripts': [
        # Main programs
        'genomes = mitty.genomes:cli',
        'reads = mitty.reads:cli',
        'god-aligner = mitty.benchmarking.god_aligner:cli',

        # Benchmarking code
        'bam2df = mitty.benchmarking.bam2df:cli',
        'evalvcf2df = mitty.benchmarking.evalvcf2df:cli',
        'edf-bdf = mitty.benchmarking.edf_bdf:cli',
        'call-diff = mitty.benchmarking.call_diff:cli',
        'df2bam = mitty.benchmarking.df2bam:cli',
        'df2vcf = mitty.benchmarking.df2vcf:cli',
        'filterdf = mitty.benchmarking.filterdf:cli',

        # Plotting functions
        'mq = mitty.benchmarking.plotting.mq:cli',

        # The stuff below may be deprecated
        'perfectbam = mitty.benchmarking.perfectbam:cli',
        'mismat = mitty.benchmarking.mismat:cli',
        'indelbam = mitty.benchmarking.indelbam:cli',
        'badbams = mitty.benchmarking.badbams:cli',
        #'alindel = mitty.benchmarking.indel_alignment_accuracy:cli',
        'indelplot = mitty.benchmarking.indel_plot:cli',
        #'aggregate-summary-report = mitty.benchmarking.aggregate_indel_plot:cli',
        'mqplot = mitty.benchmarking.mq_plot:cli',
        'dplot = mitty.benchmarking.d_plot:cli',
        'aligner-summary-report = mitty.benchmarking.aligner_analysis_summary:cli',
        'vcf2pop = mitty.lib.vcf2pop:cli',
        'bam2tfq = mitty.benchmarking.convert_bam_to_truth_fastq:cli',
        #'alindel_plot = mitty.benchmarking.indel_alignment_accuracy_plot:cli',
        'misplot = mitty.benchmarking.misalignment_plot:cli',
        'aggregate-data = mitty.benchmarking.aggregate_data:cli',
        #'acubam = mitty.benchmarking.bam_accuracy:cli',
        'migratedb = mitty.util.db_migrate:cli',
        'plot_gc_bias = mitty.util.plot_gc_bias:cli',
        'splitta = mitty.util.splitta:cli',
        'kmers = mitty.util.kmers:cli',
        'pybwa = mitty.util.pybwa:cli']
    },
    install_requires=[
      'cython',
      'setuptools>=11.0.0',
      'numpy>=1.9.0',
      'docopt>=0.6.2',
      'click>=3.3',
      'pysam==0.8.4',  # 0.9.0 gives a StringIO error. 0.8.3 is confused by some set_tag operations
      'h5py>=2.5.0',
      'matplotlib>=1.3.0',
      'scipy',
      'pandas>=0.18.1',
      'tables>=3.2.2'
    ],
)