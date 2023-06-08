import argparse, sys, StringIO
import os
import numpy as np
import pysam

ar=os.path.dirname(os.path.realpath(__file__)).split('/')
svtpath='/'.join(ar[0:(len(ar)-1)])
sys.path.insert(1, svtpath)
from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant
from collections import namedtuple
import svtools.utils as su


def add_arguments_to_parser(parser):
    parser.add_argument('-i', '--vcf', metavar='<VCF>', dest='manta_vcf', help="manta input vcf")
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true')

def command_parser():
    parser = argparse.ArgumentParser(description="fix up manta vcfs for lsort and lmerge")
    add_arguments_to_parser(parser)
    return parser

def convert_dup(v):
    v.info['SVTYPE']='INS'
    v.alt='<INS>';
    return v
        
def run_from_args(args):

  vcf = Vcf()
  vcf_out=sys.stdout
  in_header = True
  header_lines = list()
  chrdict = {}
  chrnum=1
  with su.InputStream(args.manta_vcf) as input_stream:
    for line in input_stream:
      if in_header:
        header_lines.append(line)
        if line[0:12] == '##contig=<ID':
          chrom=line.replace(',', '=').split('=')[2]
          chrdict[chrom]=chrnum
          chrnum=chrnum+1
        if line[0:6] == '#CHROM':
          in_header=False
          vcf.add_header(header_lines)
          vcf.add_info('NCALLERS', '1', 'Integer', 'Number of callers')
          vcf.add_info('CALLERS', '1', 'String', 'Callers')
          vcf.add_info('HITS_STR', '1', 'Integer', 'binary hits str')
          vcf.add_info('HITS_SEGDUP', '1', 'Integer', 'binary hits segdup')
          vcf.add_info('HITS_SAT', '1', 'Integer', 'binary hits sat')
          vcf.add_info('HAS_ILL', '1', 'Integer', 'binary has illumina support')
          vcf.add_info('HAS_LR', '1', 'Integer', 'binary has lr support')
          vcf.add_info('HAS_ASS', '1', 'Integer', 'binary has assembly support')
          vcf.add_info('HAS_MAT', '1', 'Integer', 'binary has support maternal')
          vcf.add_info('HAS_PAT', '1', 'Integer', 'binary has support paternal')
          vcf.add_info('HAS_DIP', '1', 'Integer', 'binary has diploid caller support')
          vcf_out.write(vcf.get_header()+'\n')
      else:
        v = Variant(line.rstrip().split('\t'), vcf)
        if v.info['SVTYPE']=='DUP':
          convert_dup(v)
        vcf_out.write(v.get_var_string()+"\n")
          

parser=command_parser()
args=parser.parse_args()
run_from_args(args)
