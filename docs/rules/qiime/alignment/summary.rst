Qiime - Phylogeny Align To Tree Mafft Fasttree
--------------------------------------------------

Generate a phylogenetic tree with align-to-tree-mafft-fasttree pipeline from the q2-phylogeny plugin.

Calls 'qiime phylogeny align-to-tree-mafft-fasttree'.

First, the pipeline uses the mafft program to perform a multiple sequence alignment of the sequences at input artifact (FeatureData[Sequence]) to create a FeatureData[AlignedSequence] QIIME 2 artifact. Next, the pipeline masks (or filters) the alignment to remove positions that are highly variable. These positions are generally considered to add noise to a resulting phylogenetic tree. Following that, the pipeline applies FastTree to generate a phylogenetic tree from the masked alignment. The FastTree program creates an unrooted tree, so in the final step in this section midpoint rooting is applied to place the root of the tree at the midpoint of the longest tip-to-tip distance in the unrooted tree.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/qiime/alignment/phylogeny_align_to_tree_mafft_fasttree.snake
- *Rule name:* qiime__phylogeny_align_to_tree_mafft_fasttree

**Input(s):**

- *rep_seqs:* representative sequences, qiime artifact of type FeatureData[Sequence], e.g. reads/qiime/joined-rep-seqs-dn-99.qza

**Output(s):**

- *alignment:* aligned sequences, qiime artifact of type FeatureData[AlignedSequence], e.g. reads/qiime/joined-aligned-rep-seqs-dn-99.qza
- *masked_alignment:* masked aligned sequences, qiime artifact of type FeatureData[AlignedSequence], e.g. reads/qiime/joined-rep-masked-aligned-seqs-dn-99.qza
- *unrooted_tree:* unrooted tree, qiime artifact of type Phylogeny[Unrooted], e.g. reads/qiime/joined-unrooted-tree-dn-99.qza
- *rooted_tree:* rooted tree, qiime artifact of type Phylogeny[Rooted], e.g. reads/qiime/joined-rooted-tree-dn-99.qza

