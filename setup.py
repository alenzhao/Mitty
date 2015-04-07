from setuptools import setup, find_packages, Extension

cy_modules = ['mitty.lib.util']


def make_extension_modules(post_processor, source_file_extension):
    return post_processor([Extension(cl, ['{:s}.{:s}'.format(cl.replace('.', '/'), source_file_extension)]) for cl in cy_modules])


def extension_modules():
    try:
        from Cython.Build import cythonize
        return make_extension_modules(cythonize, 'pyx')
    except ImportError:
        return make_extension_modules(lambda arg: arg, 'c')


setup(
    name='mitty',
    version='1.1.0.dev',
    description='Simulator for genomic data',
    author='Seven Bridges Genomics',
    author_email='kaushik.ghose@sbgenomics.com',
    packages=find_packages(include=['mitty*']),
    include_package_data=True,
    entry_points={
      # Register the built in plugins
      'mitty.plugins.sfs': ['double_exp = mitty.plugins.site_frequency.double_exp'],
      'mitty.plugins.variants': ['snp = mitty.plugins.variants.snp_plugin',
                                 #'delete = mitty.plugins.variants.delete_plugin',
                                 #'bounded_delete = mitty.plugins.variants.bounded_len_delete_plugin',
                                 #'insert = mitty.plugins.variants.insert_plugin',
                                 #'inversion = mitty.plugins.variants.inversion_plugin',
                                 #'low_entropy_insert = mitty.plugins.variants.low_entropy_insert_plugin'
                                 ],
      'mitty.plugins.reads': ['simple_sequential = mitty.plugins.reads.simple_sequential_plugin',
                              'simple_illumina = mitty.plugins.reads.simple_illumina_plugin'],
      # Register example tool wrapper
      'mitty.benchmarking.tools': ['bwa = mitty.benchmarking.tool_wrappers.bwa'],
      # Command line scripts
      'console_scripts': ['genomes = mitty.genomes:cli',
                          'reads = mitty.reads:cli',
                          'inspect = mitty.util.inspect:cli',
                          'perfectbam = mitty.util.perfectbam:cli',
                          'splitta = mitty.util.splitta:cli',
                          'checkbam = mitty.util.checkbam:cli']
    },
    install_requires=[
      'setuptools>=0.7',
      'numpy>=1.9.0',
      'docopt>=0.6.2',
      'pysam>=0.8.1',
      'PyVCF==0.7.0dev'
    ],
    ext_modules=extension_modules(),
)
