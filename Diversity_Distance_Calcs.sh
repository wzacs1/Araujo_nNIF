#!/bin/bash

# Requires: QIIME2 version 2021.8 and the following input files processed from raw sequences:

# Starts with QIIME2 processed files:
#    qza/repseq.qza
#    qza/table_bySeqID.qza
#    qza/taxonomy.qza
#    qza/tree_root.qza
# And a metadata ("map") file relating sample IDs to groupings:
#    sampleSheet.txt

# Make some directories and set variables for map file and number processes
mkdir -p tax
mkdir -p div
mkdir -p q2_viz
MAP=$PWD/sampleSheet.txt
NJOBS=4

# Remove ASVs that do not have a phylum call at least. We can tell these are non-16S mouse off-target sequences, enriched in peritoneal fluid samples.
qiime taxa filter-table --i-table qza/table_bySeqID.qza --o-filtered-table qza/table_bySeqID_wPhyla.qza --i-taxonomy qza/taxonomy.qza --p-include p__
qiime feature-table summarize --i-table qza/table_bySeqID_wPhyla.qza --o-visualization q2_viz/table_bySeqID_wPhyla.qzv --m-sample-metadata-file $MAP
qiime feature-table filter-seqs --i-data qza/repseq.qza --i-table qza/table_bySeqID_wPhyla.qza --o-filtered-data qza/repseq_wPhyla.qza
qiime phylogeny filter-tree --i-table qza/table_bySeqID_wPhyla.qza --i-tree qza/tree_root.qza --o-filtered-tree qza/tree_root_wPhyla.qza

# Create an ASV table with just stool samples and an ASV table with just peritoneal flood ("PF"):
qiime feature-table filter-samples --i-table qza/table_bySeqID_wPhyla.qza --o-filtered-table qza/table_bySeqID_wPhyla_PF.qza --m-metadata-file $MAP --p-where "[Organ] == 'PF'"
qiime feature-table summarize --i-table qza/table_bySeqID_wPhyla_PF.qza --m-sample-metadata-file $MAP --o-visualization q2_viz/table_bySeqID_wPhyla_PF.qzv
qiime feature-table filter-samples --i-table qza/table_bySeqID_wPhyla.qza --o-filtered-table qza/table_bySeqID_wPhyla_Stool.qza --m-metadata-file $MAP --p-where "[Organ] == 'Stool'"
qiime feature-table summarize --i-table qza/table_bySeqID_wPhyla_Stool.qza --m-sample-metadata-file $MAP --o-visualization q2_viz/table_bySeqID_wPhyla_Stool.qzv

# Run the core phylogenetic diversity pipeline on the 3 different ASV tables (All, PF, Stool), passing minimum # observations in each group:
# Figure 4C from PCoA created with "All" sample set
qiime diversity core-metrics-phylogenetic --i-table qza/table_bySeqID_wPhyla_PF.qza --i-phylogeny qza/tree_root.qza --output-dir div/PF --m-metadata-file $MAP --p-n-jobs-or-threads $NJOBS --p-sampling-depth 4476
qiime diversity core-metrics-phylogenetic --i-table qza/table_bySeqID_wPhyla_Stool.qza --i-phylogeny qza/tree_root.qza --output-dir div/Stool --m-metadata-file $MAP --p-n-jobs-or-threads $NJOBS --p-sampling-depth 160198
qiime diversity core-metrics-phylogenetic --i-table qza/table_bySeqID_wPhyla.qza --output-dir div/All --p-sampling-depth 4476 --i-phylogeny qza/tree_root.qza --m-metadata-file $MAP --p-n-jobs-or-threads $NJOBS


# Test if treatment groups are significantly different (Figure 4D):
for Organ in div/Stool div/PF
do

  for DivMetric in bray_curtis jaccard unweighted_unifrac weighted_unifrac
  do
    qiime diversity beta-group-significance --i-distance-matrix ${Organ}/${DivMetric}_distance_matrix.qza --m-metadata-file $MAP --m-metadata-column Treatment --p-permutations 9999 --p-pairwise --o-visualization ${Organ}/${DivMetric}_TreatSig.qzv
  done

done

# We want to know the dissimilarity of stool and peritoneal fluid in each treatment group. To quickly calculate these (exported from resultant .qzv files) and ensure the same ASV table is used (so no differences due to different rarefactions),
#   subset the rarefied table from core diversity calculations to each treatment group and then run beta group significance command (calculates distances in easily exportable format which is then replotted in Prism software).
for treat in NIF NIF_Meropenem Vehicle Meropenem
do
  qiime feature-table filter-samples --i-table div/All/rarefied_table.qza --o-filtered-table div/All/rarefied_table_${treat}.qza --m-metadata-file $MAP --p-where "[Treatment] == '${treat}'"
done

for table in div/All/rarefied_table_NIF.qza div/All/rarefied_table_NIF_Meropenem.qza div/All/rarefied_table_Vehicle.qza div/All/rarefied_table_Meropenem.qza
do
  qiime diversity core-metrics-phylogenetic --i-table ${table} --i-phylogeny qza/tree_root.qza --m-metadata-file $MAP --p-sampling-depth 4476 --output-dir ${table%.qza}_div --p-n-jobs-or-threads $NJOBS
done

for treatgroup in div/All/rarefied_table_NIF_div div/All/rarefied_table_NIF_Meropenem_div div/All/rarefied_table_Vehicle_div div/All/rarefied_table_Meropenem_div
do
  qiime diversity beta-group-significance --i-distance-matrix ${treatgroup}/bray_curtis_distance_matrix.qza --m-metadata-file $MAP --m-metadata-column Organ --p-permutations 9999 --o-visualization ${treatgroup}/bray_curtis_distance_OrganSig
done
