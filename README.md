# Files Descriptions:

- [sampleSheet.txt](sampleSheet.txt): Sample map relating sample metadata to sample IDs.  
- [dataProcessingScripts/](dataProcessingScripts/): Directory containing nextflow and bash scripts to process and quality check fastq sequence files and produce QIIME2 artefact files. 
- [qza/](qza/): Directory containing QIIME2 artefact files (repseq.qza, tree_root.qza, taxonomy.qza and table_bySeqID.qza) resulting from initial processing of fastq sequence files and used for further analyses.  
- [Diversity_Distance_Calcs.sh](Diversity_Distance_Calcs.sh): Bash script using QIIME2 commands to calculate diversity stats including distance values, PCoA and PERMANOVA significance values from Figure 4C-E.   
- [ANCOM-BC_markdown.Rmd](ANCOM-BC_markdown.Rmd): R markdown with commands used for running ANCOM-BC for differential taxa abundance.  
