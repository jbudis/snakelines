import gzip
import sys

in_fastq = sys.argv[1]
out_txt = sys.argv[2]
sid = sys.argv[3]

readgroup = ''

with gzip.open(in_fastq, 'rt') as fastq_file:
    header = fastq_file.readline()
    
    if header is not None:
        # assume illumina header
        # @AP2-11:127:H53WFDMXX:1:1101:1271:1031 1:N:0:GGACTCCT+GTAAGGAG
        segments = header.split(':')
        assume_illumina = True
        if len(segments) >= 3:
            flowcell = segments[2]
            platform = 'ILLUMINA'
        else:
            # assume bgi header
            # @ERR2618717.1 CL200036657L2C001R002_104504 length=150
            segments = header.split(' ')
            if len(segments) >= 2:
                flowcell = segments[1]
                platform = 'BGI'
            else:
                raise ValueError('unable to obtain flowcell from read header: ' % header)
        
        readgroup += '@RG\t'
        readgroup += 'ID:%s.%s\t' % (flowcell, sid)
        readgroup += 'LB:%s.%s\t' % (flowcell, sid)
        readgroup += 'PL:%s\t' % platform
        readgroup += 'SM:%s' % sid
    else:
        raise ValueError('input FASTQ %s is empty' % header)

if not readgroup:
    raise ValueError('readgroup is empty')

with open(out_txt, 'wt') as out_file:
    out_file.write(readgroup + '\n')
