import sys
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio import SeqIO

input_alignment, input_tree = sys.argv[1:3]
threads = sys.argv[4]
output_svg, output_ascii = sys.argv[5:7]

# If all sequences have defined proportion of reads they take in sample, append it to OTU labels
labels = {}
for contig in SeqIO.parse(input_alignment, 'fasta'):
    ratio_items = re.findall(r"ratio=(0\.\d+|1)", contig.description)
    if ratio_items and len(ratio_items) == 1:
        ratio = float(ratio_items[0])
        labels[contig.id] = '{contig.id}: {ratio: .2f}%'.format(contig=contig, ratio=ratio*100)
    else:
        labels = None
        break

label_func = (lambda clade: labels.get(clade.name, '')) if labels else str

phylo_tree = Phylo.read(input_tree, "newick")
phylo_tree.ladderize()
with open(output_ascii, 'w') as out_ascii:
    Phylo.draw_ascii(phylo_tree, file=out_ascii)
Phylo.draw(phylo_tree, label_func=label_func, do_show=False)
plt.savefig(output_svg,format='svg', bbox_inches='tight', dpi=300)

