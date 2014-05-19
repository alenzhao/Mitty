"""This is the wrapper for vcf2seq.py tool.
Command line parameters are

<ref_seq>  <var_seq> <chrom> <vcf_file> [--ploidy=PL] [--block_size=BS] [-v]

"""
from sbgsdk import define, Process, require
import os


@require(mem_mb=4048, cpu=require.CPU_SINGLE)
class Vcf2Seq(define.Wrapper):
  class Inputs(define.Inputs):
    ref_seq = define.input(name='Reference', description='The reference sequence (.smalla format)', required=True)
    vcf_file = define.input(name='VCF file', description='VCF file', required=True)

  class Outputs(define.Outputs):
    var_seqs = define.output(name='Output sequence (s)',
      description='smalla file(s) containing variant sequence(s) with listed variants', list=True)

  class Params(define.Params):
    var_seq_prefix = define.string(required=True, default='variant_seq',
                                   description='Prefix for the ouput smalla file(s)')
    chromosome = define.string(required=True)
    ploidy = define.integer(required=True, default=1, min=1, max=2,
                            description='How many sets of chromosomes do we have? We can take a diploid VCF file'
                                        '(with genotype info) and discard the genotype information - treat all'
                                        'mutations as homozygous - by putting in ploidy=1 here')

  def execute(self):
    #output_name = self.params.output_vcf_name or input_name + '_variants.vcf'
    #output_name = input_name + '_variants.vcf'
    output_dir = 'OUTPUT'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)

    p = Process('python', '/Mitty/vcf2seq.py',
                self.inputs.ref_seq,
                os.path.join(output_dir, self.params.var_seq_prefix),
                self.params.chromosome,
                self.inputs.vcf_file,
                '--ploidy', self.params.ploidy, '-v')
    p.run()

    self.outputs.var_seqs = [os.path.join(output_dir,
                                          self.params.var_seq_prefix + '_' + str(n) + '.smalla') for n in range(self.params.ploidy)]
    self.outputs.var_seqs.meta = self.outputs.var_seqs.make_metadata(file_type='smalla', type='list')


def test_vcf2seq():
  """Test with simple VCF file and the porcine circovirus test data"""
  from vcf2seq_wrapper import Vcf2Seq
  import pysam  # Need this to compress and index our sample VCF file

  # Create a test vcf file. You may recognize this from the Readme
  with open('/sbgenomics/test-data/variants.vcf', 'w') as f:
    f.write("""##fileformat=VCFv4.1
##fileDate=...
##source=...
##reference=Porcine circovirus
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample
1\t99\t.\tC\tG\t96\tPASS\t.\tGT\t1/1
1\t102\t.\tCCGTTACCGC\tC\t96\tPASS\t.\tGT\t1/1
1\t187\t.\tA\tC\t96\tPASS\t.\tGT\t1/1
1\t196\t.\tTCCTGGG T\t96\tPASS\t.\tGT\t1/1
1\t277\t.\tT\tC\t96\tPASS\t.\tGT\t1/1
1\t283\t.\tTACTACAG\tT\t96\tPASS\t.\tGT\t1/1
1\t361\t.\tAGTGCTGTTA\tA\t96\tPASS\t.\tGT\t1/1
1\t374\t.\tT\tA\t96\tPASS\t.\tGT\t1/1
1\t465\t.\tCTACCACTCC\tC\t96\tPASS\t.\tGT\t1/1
1\t570\t.\tT\tC\t96\tPASS\t.\tGT\t1/1
1\t657\t.\tA\tG\t96\tPASS\t.\tGT\t1/1
1\t663\t.\tCAGAGAA C\t96\tPASS\t.\tGT\t1/1
""")
  pysam.tabix_compress('/sbgenomics/test-data/variants.vcf', '/sbgenomics/test-data/variants.vcf.gz', force=True)
  pysam.tabix_index('/sbgenomics/test-data/variants.vcf.gz', force=True, preset='vcf')

  # This is the mutated sequence
  expected_mut_seq = 'ATGACGTATCCAAGGAGGCGTTACCGGAGAAGAAGACACCGCCCCCGCAGCCATCTTGGCCAGATCCTCCGCCGCCGCCCCTGGCTCGTCCACCCCCGGCACTGGAGAAGGAAAAACGGCATCTTCAACACCCGCCTCTCCCGCACCTTCGGATATACTATCAAGCGAACCACAGTCCAAACGCCCTCGGTGGACATGATGAGATTCAATATTAATGACTTTCTTCCCCCAGGAGGGGGCTCAAACCCCCGCTCTGTGCCCCTTGAATAATAAGAAAGGTTAAGGTTGAATTCTGGCCCTGCTCCCCGATCACCCAGGGTGACAGGGGAGTGGGCTCCATTCAAGATGATAACTTTGTAACAAAGGCCACAGCCCTCACCTATGACCCCTATGTAAACTACTCCTCCCGCCATACCATAACCCAGCCCTTCTCCCGCTACTTTACCCCCAAACCTGTCCTAGATTCCACTATTGATTACTTCCAACCAAACAACAAAAGAAATCAGCTGTGGCTGAGACTACAAACTGCCGGAAATGTAGACCACGTAGGCCTCGGCACTGCGTTCGAAAACAGTATATACGACCAGGAATACAATATCCGTGTAACCATGTATGTGCAATTCTTTAATCTTAAAGACCCCCCACTTAACCCTTAG'

  inputs = {'ref_seq': '/sbgenomics/test-data/porcine_circovirus.smalla',
            'vcf_file': '/sbgenomics/test-data/variants.vcf'}
  params = {'var_seq_prefix': 'variant_seq',
            'chromosome': '1',
            'ploidy': 2}
  wrp = Vcf2Seq(inputs, params)
  outputs = wrp.test()

  # Test to see the written sequences are correct
  for fname in outputs.var_seqs:
    with open(fname, 'r') as f:
      assert f.read() == expected_mut_seq
