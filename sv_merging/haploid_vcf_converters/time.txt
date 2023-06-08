import argparse, sys, StringIO
import os
#import numpy as np
#import pysam

from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant
from collections import namedtuple
import svtools.utils as su


def add_arguments_to_parser(parser):
    parser.add_argument('-i', '--vcf', metavar='<VCF>', dest='manta_vcf', help="manta input vcf")
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true')
    parser.add_argument('-s', '--slop', dest='slop',  default=0,  required=False, help='padding to either side')
    parser.add_argument('-c', '--caller', dest='caller', required=True, help='variant caller')
    parser.add_argument('-h', '--haplotype', dest='haplotype', required=True)


chrdict={
"chr1":1,
"chr2":2,
"chr3":3,
"chr4":4,
"chr5":5,
"chr6":6,
"chr7":7,
"chr8":8,
"chr9":9,
"chr10":10,
"chr11":11,
"chr12":12,
"chr13":13,
"chr14":14,
"chr15":15,
"chr16":16,
"chr17":17,
"chr18":18,
"chr19":19,
"chr20":20,
"chr21":21,
"chr22":22,
"chrX":23,
"chrY":24,
"chrM":25,
"chr1_KI270706v1_random":26,
"chr1_KI270707v1_random":27,
"chr1_KI270708v1_random":28,
"chr1_KI270709v1_random":29,
"chr1_KI270710v1_random":30,
"chr1_KI270711v1_random":31,
"chr1_KI270712v1_random":32,
"chr1_KI270713v1_random":33,
"chr1_KI270714v1_random":34,
"chr2_KI270715v1_random":35,
"chr2_KI270716v1_random":36,
"chr3_GL000221v1_random":37,
"chr4_GL000008v2_random":38,
"chr5_GL000208v1_random":39,
"chr9_KI270717v1_random":40,
"chr9_KI270718v1_random":41,
"chr9_KI270719v1_random":42,
"chr9_KI270720v1_random":43,
"chr11_KI270721v1_random":44,
"chr14_GL000009v2_random":45,
"chr14_GL000225v1_random":46,
"chr14_KI270722v1_random":47,
"chr14_GL000194v1_random":48,
"chr14_KI270723v1_random":49,
"chr14_KI270724v1_random":50,
"chr14_KI270725v1_random":51,
"chr14_KI270726v1_random":52,
"chr15_KI270727v1_random":53,
"chr16_KI270728v1_random":54,
"chr17_GL000205v2_random":55,
"chr17_KI270729v1_random":56,
"chr17_KI270730v1_random":57,
"chr22_KI270731v1_random":58,
"chr22_KI270732v1_random":59,
"chr22_KI270733v1_random":60,
"chr22_KI270734v1_random":61,
"chr22_KI270735v1_random":62,
"chr22_KI270736v1_random":63,
"chr22_KI270737v1_random":64,
"chr22_KI270738v1_random":65,
"chr22_KI270739v1_random":66,
"chrY_KI270740v1_random":67,
"chrUn_KI270302v1":68,
"chrUn_KI270304v1":69,
"chrUn_KI270303v1":70,
"chrUn_KI270305v1":71,
"chrUn_KI270322v1":72,
"chrUn_KI270320v1":73,
"chrUn_KI270310v1":74,
"chrUn_KI270316v1":75,
"chrUn_KI270315v1":76,
"chrUn_KI270312v1":77,
"chrUn_KI270311v1":78,
"chrUn_KI270317v1":79,
"chrUn_KI270412v1":80,
"chrUn_KI270411v1":81,
"chrUn_KI270414v1":82,
"chrUn_KI270419v1":83,
"chrUn_KI270418v1":84,
"chrUn_KI270420v1":85,
"chrUn_KI270424v1":86,
"chrUn_KI270417v1":87,
"chrUn_KI270422v1":88,
"chrUn_KI270423v1":89,
"chrUn_KI270425v1":90,
"chrUn_KI270429v1":91,
"chrUn_KI270442v1":92,
"chrUn_KI270466v1":93,
"chrUn_KI270465v1":94,
"chrUn_KI270467v1":95,
"chrUn_KI270435v1":96,
"chrUn_KI270438v1":97,
"chrUn_KI270468v1":98,
"chrUn_KI270510v1":99,
"chrUn_KI270509v1":100,
"chrUn_KI270518v1":101,
"chrUn_KI270508v1":102,
"chrUn_KI270516v1":103,
"chrUn_KI270512v1":104,
"chrUn_KI270519v1":105,
"chrUn_KI270522v1":106,
"chrUn_KI270511v1":107,
"chrUn_KI270515v1":108,
"chrUn_KI270507v1":109,
"chrUn_KI270517v1":110,
"chrUn_KI270529v1":111,
"chrUn_KI270528v1":112,
"chrUn_KI270530v1":113,
"chrUn_KI270539v1":114,
"chrUn_KI270538v1":115,
"chrUn_KI270544v1":116,
"chrUn_KI270548v1":117,
"chrUn_KI270583v1":118,
"chrUn_KI270587v1":119,
"chrUn_KI270580v1":120,
"chrUn_KI270581v1":121,
"chrUn_KI270579v1":122,
"chrUn_KI270589v1":123,
"chrUn_KI270590v1":124,
"chrUn_KI270584v1":125,
"chrUn_KI270582v1":126,
"chrUn_KI270588v1":127,
"chrUn_KI270593v1":128,
"chrUn_KI270591v1":129,
"chrUn_KI270330v1":130,
"chrUn_KI270329v1":131,
"chrUn_KI270334v1":132,
"chrUn_KI270333v1":133,
"chrUn_KI270335v1":134,
"chrUn_KI270338v1":135,
"chrUn_KI270340v1":136,
"chrUn_KI270336v1":137,
"chrUn_KI270337v1":138,
"chrUn_KI270363v1":139,
"chrUn_KI270364v1":140,
"chrUn_KI270362v1":141,
"chrUn_KI270366v1":142,
"chrUn_KI270378v1":143,
"chrUn_KI270379v1":144,
"chrUn_KI270389v1":145,
"chrUn_KI270390v1":146,
"chrUn_KI270387v1":147,
"chrUn_KI270395v1":148,
"chrUn_KI270396v1":149,
"chrUn_KI270388v1":150,
"chrUn_KI270394v1":151,
"chrUn_KI270386v1":152,
"chrUn_KI270391v1":153,
"chrUn_KI270383v1":154,
"chrUn_KI270393v1":155,
"chrUn_KI270384v1":156,
"chrUn_KI270392v1":157,
"chrUn_KI270381v1":158,
"chrUn_KI270385v1":159,
"chrUn_KI270382v1":160,
"chrUn_KI270376v1":161,
"chrUn_KI270374v1":162,
"chrUn_KI270372v1":163,
"chrUn_KI270373v1":164,
"chrUn_KI270375v1":165,
"chrUn_KI270371v1":166,
"chrUn_KI270448v1":167,
"chrUn_KI270521v1":168,
"chrUn_GL000195v1":169,
"chrUn_GL000219v1":170,
"chrUn_GL000220v1":171,
"chrUn_GL000224v1":172,
"chrUn_KI270741v1":173,
"chrUn_GL000226v1":174,
"chrUn_GL000213v1":175,
"chrUn_KI270743v1":176,
"chrUn_KI270744v1":177,
"chrUn_KI270745v1":178,
"chrUn_KI270746v1":179,
"chrUn_KI270747v1":180,
"chrUn_KI270748v1":181,
"chrUn_KI270749v1":182,
"chrUn_KI270750v1":183,
"chrUn_KI270751v1":184,
"chrUn_KI270752v1":185,
"chrUn_KI270753v1":186,
"chrUn_KI270754v1":187,
"chrUn_KI270755v1":188,
"chrUn_KI270756v1":189,
"chrUn_KI270757v1":190,
"chrUn_GL000214v1":191,
"chrUn_KI270742v1":192,
"chrUn_GL000216v2":193,
"chrUn_GL000218v1":194,
"chrEBV":195}

def command_parser():
    parser = argparse.ArgumentParser(description="fix up vcfs for lsort and lmerge")
    add_arguments_to_parser(parser)
    return parser


def convert_variant(v,  chrdict):
     set_read_counts(v)
     set_cis_prs(v)
     if v.get_info('SVTYPE')=='DEL':
         convert_del(v)
     elif v.get_info('SVTYPE')=='DUP:INT' or v.get_info('SVTYPE')=='DUP:TANDEM':
         convert_dup(v)
     elif v.get_info('SVTYPE')=='INV':
         convert_inv(v)
     elif v.get_info('SVTYPE')=='INS':
         convert_ins(v)
     elif v.get_info('SVTYPE')=='BND':
         convert_bnd(v, chrdict)
        
def split_ci(ci):
     return[int(ci.split(',')[0]),  int(ci.split(',')[1])]

def uniform_pr(length):
    pr=np.ones(length, dtype='float64')/length
    pr1=','.join( map(str, pr))
    return pr1

def set_read_counts(var):
    pe=2
    sr=2
    var.info['PE']=pe
    var.info['SR']=sr
    var.info['SU']=pe+sr

def set_cis_prs(v):
    imprec=True
    cipos='0,0'; ciend='0,0'
    prpos=1.0; prend=1.0
    v.info['CIPOS']=cipos
    v.info['CIEND']=ciend
    v.info['CIPOS95']=cipos
    v.info['CIEND95']=ciend
    v.info['PRPOS']=prpos
    v.info['PREND']=prend
    v.set_info('IMPRECISE', imprec)
    
def convert_del(var):
    var.alt='<DEL>'
    var.info['STRANDS']='+-:'+str(var.info['SU'])
    var.ref='N'

def convert_dup(var):
    var.ref='N'
    var.alt='<DUP>'
    var.info['SVTYPE']='DUP'
    var.info['STRANDS']='-+:'+str(var.info['SU'])
    var.info['END']=var.pos+var.get_info('SVLEN')
    var.info['SVTYPE_ORIG']='DUP'


def convert_inv(var):
    var.ref='N'
    var.alt='<INV>'
    strands='++:2,--:2'
    var.info['STRANDS']=strands+str(var.info['SU'])
    var.info['SVLEN']=int(var.info['END'])-var.pos

def convert_ins(var):
    var.ref='N'
    var.alt='<DUP>'
    var.info['STRANDS']='-+:'+str(var.info['SU'])
    var.info['SVTYPE']='DUP'
    var.info['SVTYPE_ORIG']='INS'
    var.info['END']=var.pos+var.get_info('SVLEN')
    
        
def convert_bnd(var, chrdict):
    var.ref='N'
    alt=var.alt
    ff=alt.find("[")
    newalt=""
    strands=""
    sep1, chrom2, breakpoint2=su.parse_bnd_alt_string(alt)

    if int(chrdict[chrom2])<int(chrdict[var.chrom]) or (chrom2==var.chrom and int(breakpoint2)<int(var.pos)):
        var.set_info('SECONDARY', True)
    if ff==0:
        strands="--:"
        ff1=alt.find("[", 1)
        newalt=alt[0:(ff1+1)]+'N'
    elif ff>0:
        strands="+-:"
        newalt='N'+alt[ff::]
    else:
        ff=alt.find("]")
        if ff==0:
            strands="-+:"
            ff1=alt.find("]", 1)
            newalt=alt[0:(ff1+1)]+'N'
        else:
            strands="++:"
            newalt='N'+alt[ff::]
    var.alt=newalt
    var.info['STRANDS']=strands+str(var.info['SU'])
        

def run_from_args(args):

  vcf = Vcf()
  vcf_out=sys.stdout
  in_header = True
  header_lines = list()
  with su.InputStream(args.manta_vcf) as input_stream:
    for line in input_stream:
      if in_header:
        header_lines.append(line)
        if line[0:6] == '#CHROM':
          in_header=False
          vcf.add_header(header_lines)
          vcf.add_info('PRPOS', '1', 'String', 'Breakpoint probability dist')
          vcf.add_info('PREND', '1', 'String', 'Breakpoint probability dist')
          vcf.add_info('STRANDS', '.', 'String', 'Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--')
          vcf.add_info('SU', '.', 'Integer', 'Number of pieces of evidence supporting the variant across all samples')
          vcf.add_info('PE', '.', 'Integer', 'Number of paired-end reads supporting the variant across all samples')
          vcf.add_info('SR', '.', 'Integer', 'Number of split reads supporting the variant across all samples')
          vcf.add_info('INSLEN_ORIG', '.', 'Integer', 'Original insertion length')
          vcf.add_info('CIPOS', '2', 'Integer', 'Confidence interval around POS for imprecise variants')
          vcf.add_info('CIEND', '2', 'Integer', 'Confidence interval around END for imprecise variants')
          vcf.add_info('CIPOS95', '2', 'Integer', 'Confidence interval (95%) around POS for imprecise variants')
          vcf.add_info('CIEND95', '2', 'Integer', 'Confidence interval (95%) around END for imprecise variants')
          vcf.add_info('SECONDARY', '0', 'Flag', 'Secondary breakend in a multi-line variant')
          vcf.add_info('SVTYPE_ORIG', '.', 'String', 'Svtype pre-conversion')
          vcf.add_info('IMPRECISE', '0', 'Flag', 'Imprecise breakpoint coordinates')
          vcf.add_info('MATEID', '.', 'String','ID of mate breakends')
          vcf.add_info('EVENT', '.', 'String', 'ID of event associated to breakend')
          vcf_out.write(vcf.get_header()+'\n')
      else:
        ll=line.rstrip().split('\t')
        ll[2]=args.caller+'_'+ll[2]
        v = Variant(ll, vcf)
        if v.get_info('SVTYPE')!='cnv':
          svlen=100
          if 'SVLEN' in v.info:
            svlen=abs(int(v.get_info('SVLEN')))
          if svlen>=50:
            convert_variant(v,  chrdict)
            vcf_out.write(v.get_var_string()+"\n")
          

parser=command_parser()
args=parser.parse_args()
run_from_args(args)
