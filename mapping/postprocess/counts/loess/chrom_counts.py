from __future__ import print_function

import argparse
import numpy as np
import os
import sys
from datetime import datetime
from subprocess import PIPE, Popen
import new_loess


def load_arguments():
    """
    Loads and parses arguments.
    :return: argparse arguments
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="chrom_counts.py -- a simple script to compute reads per chromosome corrected by loess GC correction")
    required = parser.add_argument_group('Required')
    required.add_argument('bam', type=nonempty_file, help='File path to SAM/BAM-file')
    required.add_argument('output', type=str, help="File path for output - file with counts.")

    parameters = parser.add_argument_group('Parameters')
    parameters.add_argument('-q', '--mapq', type=positive_nonzero_int, default=1, help="Minimal map quality. Default=1")
    parameters.add_argument('--flen-min', type=positive_nonzero_int, default=1, help="Minimal fragment length. Default=1")
    parameters.add_argument('--flen-max', type=positive_nonzero_int, default=400, help="Maximal fragment length. Default=400")
    parameters.add_argument('--hg38', action="store_true", help="Use hg38 genome instead hg19.")
    parameters.add_argument('-l', '--loess', action="store_true", help="Do loess correction.")

    args = parser.parse_args()

    return args


def positive_int(value, max_limit=None):
    """
    Represents positive decimal number, 0 included
    :param value: string value to estimate
    :param max_limit: maximal allowed value, skip validation if None
    :return: integer value, if can be converted into positive int else ArgumentTypeError
    """
    try:
        int_value = int(value)
    except ValueError:
        error = "Value %s is not integer" % value
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    if int_value < 0:
        error = "Value %s is not positive integer" % value
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    if max_limit and max_limit < int_value:
        error = "Value %s must be lower or equal to %s" % (value, max_limit)
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    return int_value


def positive_nonzero_int(value):
    """
    Represents positive decimal number, 0 excluded
    :param value: string value to estimate
    :return: integer value, if can be converted into positive int else ArgumentTypeError
    """
    int_value = positive_int(value)
    if int_value == 0:
        error = "Value %s cannot be 0" % value
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    return int_value


def nonempty_file(file_path):
    """
    Checks if the filename in input is non-empty file.
    :param file_path: str - filename to check
    :return: str - filename
    """
    if not os.path.exists(file_path):
        error = "File %s does not exist" % file_path
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    if os.path.getsize(file_path) == 0:
        error = "File %s is empty" % file_path
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    return file_path


CHROM_NAMES = {"chr%d" % (n + 1): n for n in range(22)}
CHROM_NAMES.update({'chrX': 22, 'chrY': 23, 'chrM': 24})
CHROM_NUMS = {v: k for k, v in CHROM_NAMES.items()}


def get_prefix(filename):
    """
    Gets a file without the suffix.
    :param filename: str - filename
    :return: str - filename without a suffix
    """
    ind = filename.rfind('.')
    if ind == -1:
        return filename
    else:
        return filename[:ind]


""" Let's go: """

# start time
start_time = datetime.now()

# read arguments
args = load_arguments()
print(sys.argv)

# data container
positions = []
chromosomes = []

# read file
if args.bam.endswith('bam'):
    to_call = 'samtools view %s' % args.bam
elif args.bam.endswith('sam'):
    to_call = 'cat %s' % args.bam
else:
    assert False, "Cannot detect the file type %s" % args.bam

p = Popen(to_call, stdout=PIPE, shell=True)
for i, line in enumerate(p.stdout):

    if line.startswith('@'):
        continue

    # progress:
    if i % 100000 == 0:
        s = "\r    processing: %8d %s, kept %8d (%.1f%%)" % (i, '{duration}'.format(duration=datetime.now() - start_time), len(positions), 0.0 if i == 0 else len(positions)/float(i)*100)
        sys.stdout.write(s)
        sys.stdout.flush()

    # count fragments (and use filters)
    split = line.strip().split()
    chr = split[2]
    pos = int(split[3])
    mapq = int(split[4])
    flen = int(split[8])

    # filter:
    if mapq < args.mapq:
        continue
    if flen < args.flen_min or flen > args.flen_max:
        continue

    # write it:
    if chr in CHROM_NAMES.keys():
        positions.append(pos)
        chromosomes.append(CHROM_NAMES[chr])

positions = np.array(positions, dtype='i4')
chromosomes = np.array(chromosomes, dtype='i2')

# save them to values.np
# values = np.array(len(positions), dtype=[('position', 'i4'), ('chromosome', 'i2')])
# values['position'] = positions
# values['chromosome'] = chromosomes
# np.save(get_prefix(args.output), values)

# apply loess
script_path = os.path.dirname(os.path.abspath(__file__))
if args.hg38:
    gc_file = '%s/gc_bins_hg38.pck' % script_path
    nn_file = '%s/nn_bins_hg38.pck' % script_path
else:
    gc_file = '%s/gc_bins_all_w.pck' % script_path
    nn_file = '%s/nn_bins_all_w.pck' % script_path

# output
if args.loess:
    loess = new_loess.get_loess_weights(chromosomes, positions, get_prefix(args.output) + '.pdf', gc_file, nn_file)
    result = {i: np.sum(loess[chromosomes == i]) for i in np.unique(chromosomes)}
else:
    result = {i: np.sum(chromosomes == i) for i in np.unique(chromosomes)}

# write result:
with open(args.output, 'w') as f:
    for k in sorted(result.keys()):
        print('%5s\t%12.3f' % (CHROM_NUMS[k], result[k]), file=f)

# write time:
end_time = datetime.now()
print('\nTime of run: {duration}'.format(duration=end_time - start_time))
