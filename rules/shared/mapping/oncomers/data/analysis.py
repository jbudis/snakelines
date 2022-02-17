import pysam
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
import numpy as np
import sys
from datetime import datetime

sid = sys.argv[1]
scriptdir = os.path.dirname(os.path.realpath(__file__))

#MELANOMA vzorky
bam = sys.argv[1]
bai = sys.argv[2]


print(f'Analysing {sid} - ' + str(datetime.now()))


if not os.path.exists(bam):
	print(f'-- Missing {bam}')
	sys.exit(0)

if not os.path.exists(bai):
	print(f'-- Missing {bai}')
	sys.exit(0)


bed = sys.argv[4]

#outdir = '/home_pfs/data/projects/population/microsat/'
#if not os.path.exists(outdir):
#	os.makedirs(outdir)

outfile = sys.argv[5]
outnpz = sys.argv[6]

if os.path.exists(outnpz):
    print(f'-- Already analysed {outnpz}')
    sys.exit(0)

samfile = pysam.AlignmentFile(bam, "rb")
num_lines = int(os.popen('wc -l ' + bed).read().split()[0])
index = 0

with open(sys.argv[3], "r") as names:
	d = {}
	j = 0
	for n in names:
		d[n.strip()] = j
		j += 1

with open(bed, "r") as ref:
	array = np.zeros(num_lines, dtype=np.uint8)	
	line = ref.readline().strip().split()
	saved = False
	start = True
	citania = iter(samfile.fetch(until_eof=True))
	while line:
		if(start or read.is_unmapped or read.mapping_quality < 2 or read.reference_end is None):
			start = False
			try:
				read = next(citania)
			except StopIteration:
				break
			else:
				continue

		if(read.reference_start <= int(line[1]) and read.reference_end >= int(line[2]) and line[0] == read.reference_name and not saved):
			beginning = False
			count = 0
			i = -1
			for x in read.get_aligned_pairs(matches_only=False, with_seq=False):
				if(x[1] is not None and int(x[1]) == int(line[1])):
					beginning = True
				if(x[0] is None):
					continue
				i += 1
				if(beginning and count > 1 and read.seq[i] != line[3]):
					if(i-count > 2 and read.qend-i > 2):
						array[index] = count
						saved = True
					break
				if(line[3] == read.seq[i]):
					count += 1
				else:
					count = 0
		if (line[0] != read.reference_name and d[line[0]] < d[read.reference_name]) or (line[0] == read.reference_name and int(line[2]) < read.reference_end):
			line = ref.readline().strip().split()
			index += 1
			saved = False
		else:
			try:
				read = next(citania)
			except StopIteration:
				break

np.savez_compressed(outfile, array)
samfile.close()

print('-- Finished - ' + str(datetime.now()))

