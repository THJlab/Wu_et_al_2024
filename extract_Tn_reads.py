import sys
import re
import pysam
from glob import glob



input_bam_file = sys.argv[1]
out_bam_file = sys.argv[2] #=input_bam_file.replace('.bam', '_Tn3x.bam')

print('extracting from ' + input_bam_file + ' to ' + out_bam_file)
def is_Tn_read(read):
    flag=read.flag
    cigar=read.cigarstring
    seq=read.seq

    if not read.is_reverse and (re.match('^..S', cigar) or re.match('^.S', cigar)): #plus strand reads (tail @start)
        tail_len=int(cigar[:cigar.index('S')])
        tail = seq[:tail_len]
        if len(tail) >= 3 and re.match('^A*$', tail):
            return True
    elif read.is_reverse and (re.search('..S$', cigar) or re.search('.S$', cigar)): #plus strand reads (tail @start)
        tail_len=int(cigar[cigar.rindex('M')+1:-1])
        tail = seq[-tail_len:]
        if len(tail) >= 3 and re.match('^T*$', tail):
            return True
    return False


in_bam = pysam.AlignmentFile(input_bam_file,"rb")
Tn_reads = (read for read in in_bam.fetch() if is_Tn_read(read))

out_bam = pysam.AlignmentFile(out_bam_file,"wb", template=in_bam)
i=0
for read in Tn_reads:
    i+=1
    x=out_bam.write(read)

print(str(i) + ' Tn tailed reads found')
in_bam.close()
out_bam.close()
