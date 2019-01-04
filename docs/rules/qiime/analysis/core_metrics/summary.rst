Qiime - Diversity Core Metrics Phylogenetic
-----------------------------------------------

Performs QIIME 2’s diversity analyses through the q2-diversity plugin.

Calls 'qiime diversity core-metrics-phylogenetic'.

Supports computing alpha and beta diversity metrics, applys related statistical tests, and generate interactive visualizations. First apply the core-metrics-phylogenetic method, which rarefies a FeatureTable[Frequency] to a user-specified depth, computes several alpha and beta diversity metrics, and generates principle coordinates analysis (PCoA) plots using Emperor for each of the beta diversity metrics. The metrics computed by default are:

* Alpha diversity
    * Shannon’s diversity index (a quantitative measure of community richness)
    * Observed OTUs (a qualitative measure of community richness)
    * Faith’s Phylogenetic Diversity (a qualitiative measure of community richness that incorporates phylogenetic relationships between the features)
    * Evenness (or Pielou’s Evenness; a measure of community evenness)

* Beta diversity
    * Jaccard distance (a qualitative measure of community dissimilarity)
    * Bray-Curtis distance (a quantitative measure of community dissimilarity)
    * unweighted UniFrac distance (a qualitative measure of community dissimilarity that incorporates phylogenetic relationships between the features)
    * weighted UniFrac distance (a quantitative measure of community dissimilarity that incorporates phylogenetic relationships between the features)

An important parameter for this script is --p-sampling-depth, which is the even sampling (i.e. rarefaction) depth. In this pipeline, sampling depth is set as Minimum of sequences count (e.g. we have 3 samples, first with 3500 sequences, 2nd with 3000 sequences and 3rd with 4500 sequences). Minimum is 3000. Value is taken from exported qzv file, exactly from reads/qiime/joined-imported/overview.html.

Note:
Because most diversity metrics are sensitive to different sampling depths across different samples, this script will randomly subsample the counts from each sample to the value provided for this parameter. For example, if you provide --p-sampling-depth 500, this step will subsample the counts in each sample without replacement so that each sample in the resulting table has a total count of 500. If the total count for any sample(s) are smaller than this value, those samples will be dropped from the diversity analysis. Choosing this value is tricky.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/qiime/analysis/core_metrics/diversity_core_metrics_phylogenetic.snake
- *Rule name:* qiime__diversity_core_metrics_phylogenetic

**Input(s):**

- *html:* html export of imported.qza file, e.g. reads/qiime/joined-imported/index.html
- *rooted_tree:* rooted tree, qiime artifact of type Phylogeny[Rooted], e.g. reads/qiime/joined-rooted-tree-dn-99.qza
- *table:* qiime artifact of type FeatureTable[Frequency], e.g. reads/qiime/joined-table.qza
- *metadata:* samples' metadata filepath in TSV format (tab-separated value), e.g. description/joined-sample-metadata.tsv

**Output(s):**

- *rarefield_table:* qiime artifact of type FeatureTable[Frequency], e.g. reads/qiime/joined-core-metrics-results/rarefied_table.qza
- *vector:* sample data for alpha diversity, qiime artifact of type SampleData[AlphaDiversity], e.g. reads/qiime/joined-core-metrics-results/shannon_vector.qza
- *distance_matrix:* distance matrix for beta diversity, qiime artifact of type DistanceMatrix, e.g. reads/qiime/joined-core-metrics-results/bray_curtis_distance_matrix.qza
- *emperor_plots:* Emperor visualization of distance matrix, e.g reads/qiime/joined-core-metrics-results/bray_curtis_emperor.qzv

