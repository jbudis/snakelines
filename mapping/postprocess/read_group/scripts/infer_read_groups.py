#!/usr/bin/python3

from argparse import ArgumentParser
from argparse import FileType

parser = ArgumentParser(description='Annotate reads with read groups infered from read ids. Script needs sample from single Illumina sequencing run')
parser.add_argument('sid', help='Sample unique identifier')
parser.add_argument('insam', type=FileType('r'), default='-')
args = parser.parse_args()

def is_header(line):
    return line.startswith('@')

def is_rg_header(line):
    return line.startswith('@RG')


TEMPLATE = '\t'.join(['@RG',
                      'ID:{rgid}', 
                      'PL:illumina', 
                      'LB:{flowcell}.{barcode}',
                      'PU:{flowcell}.{lane}',
                      'SM:{sid}'])

def infer_read_groups(line):

    rid = line[:line.find('\t')]
    items = rid.split(':')
    flowcell = items[2]
    prefix = ':'.join(items[:3])
    
    rgids = {}
    rgs = []
    for lane in [1,2,3,4]:
        rgid = '%s.l%d' % (args.sid, lane)
        rg = TEMPLATE.format(rgid=rgid, sid=args.sid, flowcell=flowcell, lane=lane, barcode=args.sid)
        rgs.append(rg)
        rgids['%s:%s' % (prefix, lane)] = rgid

    return rgids, rgs, prefix

is_first_read = True

for line in args.insam:

    # remove original read groups
    if is_rg_header(line):
        pass
    
    # keep all other headers unchanged
    elif is_header(line):
        print(line, end='')

    # all headers written, infer and write read group before first read
    else:
        if is_first_read:
    
            rgids, rgs, prefix = infer_read_groups(line)

            # print @RG headers
            for rg in rgs:
                print(rg)

            is_first_read = False
        
        rgid = rgids[line[:len(prefix) + 2]]
        last_tag_start = line.rfind('\t')
        has_rg_tag = line[last_tag_start+1] == 'R' and line[last_tag_start+2] == 'G'
        cut_pos = last_tag_start if has_rg_tag else -1
        
        print(line[:cut_pos], end='\t')
        print('RG:Z:', end='')
        print(rgid)
