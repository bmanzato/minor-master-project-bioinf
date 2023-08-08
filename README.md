# Comparing prioritized cell types from GWAS using single nucleus-RNA or -ATAC sequencing (in Alzheimer individuals and healthy controls)
## Minor Research Project for the Master in Bioinformatics and System Biology (VU Amsterdam)


The main aim of the project is develop a pipeline to link scRNAseq and scATACseq data to GWAS summary statistics.

The input are two count matrices: Gene Expression matrix for scRNAseq and Peak intensity matrix for scATACseq from [Morabito et al.2021](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174367) Features in the scRNAseq datasets are genes, in scATACseq are cis-candidate regulatory elements (cCREs). 

```preprocessing_rna.py``` and ```preprocessing_atac.py``` filter the intial datasets, divides Cases (AD) and Controls and averages the count matrix by grouping cells belonging to the same cell type.

Next, we used the [EWCE package](https://github.com/NathanSkene/EWCE) to create a specificity matrix (one for each sc technology dataset) for each feature in each cell type. Script: ```specificity.R```.

```specificity_analysis.py``` creates the Feature Sets for each cell type. In this scripts two methods are utilized to egnerate the Sets:
1. Top 10% most specific features: after filtering out features with zero value as specificity, we took the top 10% of most specific features for each cell type.
2. One-cell-type specific feature: we plotted the distribution of specificity scores for each feature to determine a threshold that could confidently identify features that are specific for one cell type only. For both scRNAseq and scATACseq features, a threshold of 0.52 was chosen to generate the set of specifically expressed/open features for each cell type.

[LDSC](https://github.com/bulik/ldsc) is used to estimate the enrichment of cell types from GWAS summary statistics (partition heritability method). ```ldsc_analysis.sh``` includes:
- creation of snp set from baseline model (v2.2)
- annotation (feature to SNP)
- compute LD scores
- run LD regression (cell-type specific)





![Flowchart](https://github.com/bmanzato/minor-master-project-bioinf/assets/74963501/3eb24462-1fcc-4c15-ad39-6dd70129b532)



