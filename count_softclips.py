import sys
import re
import pysam
from glob import glob
from collections import Counter


input_bam_file = sys.argv[1]
out_file = sys.argv[2] #=input_bam_file.replace('.bam', '_Tn3x.bam')

print('extracting from ' + input_bam_file + ' to ' + out_file)

RC={'A':'T','C':'G','G':'C','T':'A','N':'N',
    'a':'t','c':'g','g':'c','t':'a','n':'n'}
def rev_com(dna):
    rc=''.join([RC[c] for c in dna[::-1]])
    return rc


def get_softclip(read):
    flag=read.flag
    cigar=read.cigarstring
    seq=read.seq
    if not read.is_reverse and (re.match('^..S', cigar) or re.match('^.S', cigar)): #plus strand reads (tail @start)
        tail_len=int(cigar[:cigar.index('S')])
        tail = rev_com(seq[:tail_len])
    elif read.is_reverse and (re.search('..S$', cigar) or re.search('.S$', cigar)): #plus strand reads (tail @start)
        tail_len=int(cigar[cigar.rindex('M')+1:-1])
        tail = seq[-tail_len:]
    else:
        tail = None
    return tail


in_bam = pysam.AlignmentFile(input_bam_file,"rb")
softclips = Counter(get_softclip(read) for read in in_bam.fetch())
in_bam.close()

zipped_cnts_list = list(zip(softclips.values(), softclips.keys()))
zipped_cnts_list=sorted(zipped_cnts_list, key=lambda x: -x[0])
with open(out_file, 'w') as outfile:
    outfile.write('\n'.join([x[1]+'\t'+str(x[0]) if x[1] is not None else 'NA\t'+str(x[0]) for x in zipped_cnts_list ]))
