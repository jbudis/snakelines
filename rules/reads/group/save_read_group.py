import gzip
import sys

in_fastq = sys.argv[1]
out_txt = sys.argv[2]
sid = sys.argv[3]

readgroup = ''

with gzip.open(in_fastq, 'rb') as fastq_file:
    line = fastq_file.readline()
    
    if line is not None:
        segments = str(line).split(':')
        if len(segments) < 3:
            raise ValueError('unable to obtain flowcell from read name')
        
        flowcell = segments[2]
        
        readgroup += '@RG\t'
        readgroup += 'ID:%s.%s\t' % (flowcell, sid)
        readgroup += 'LB:%s.%s\t' % (flowcell, sid)
        readgroup += 'PL:ILLUMINA\t'
        readgroup += 'SM:%s' % sid
    else:
        raise ValueError('input FASTQ %s is empty' % line)

if not readgroup:
    raise ValueError('readgroup is empty')

with open(out_txt, 'wt') as out_file:
    out_file.write(readgroup + '\n')
