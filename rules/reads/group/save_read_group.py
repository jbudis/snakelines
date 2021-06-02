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
        if len(segments) >= 3:
            flowcell = segments[2]
            platform = 'ILLUMINA'
        else:
            # assume bgi header
            segments = header.split(' ')
            if len(segments) >= 2:
                # stupid BGI header ?
                # @ERR2618717.1 CL200036657L2C001R002_104504 length=150
                bgi_header = segments[1]
            else:
                # genuine BGI header, remove @
                bgi_header = segments[0][1:]

            # upper case is standard
            bgi_header = bgi_header.upper()

            if bgi_header[0] == 'C':
                # CL200036657L2C001R002_104504
                platform = 'BGISEQ-500'
                flowcell = bgi_header[:13]
            elif bgi_header[0] == 'V':
                # V300038198L4C001R0010019425/1
                platform = 'MGISEQ-2000'  # aka DNBSEQ-G400
                flowcell = bgi_header[:12]
            else:
                platform = 'unknown'
                flowcell = 'unknown'
                # raise ValueError('unknown format of read header: %s' % header)


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
