# Cholangiocarcinoma: Asian vs Western Differential Expression

Comparison of tumor vs normal gene expression in cholangiocarcinoma
(a rare bile duct cancer) across two populations, to find genes that
are consistently changed in both.

## Data
- Asian: GSE107943 (intrahepatic CCA), 18 paired tumor/normal samples, raw counts
- Western: TCGA-CHOL (UCSC Xena), 9 paired tumor/normal samples, log2 normalized

## Methods
- Asian: DESeq2 (paired design: ~ patient + condition)
- Western: limma (paired design)
- Thresholds: adj. p < 0.05, |log2FC| > 1
- Overlap of the two DEG lists, filtered for the same direction of change

## Results
- Asian DEGs: 10,267 | Western DEGs: 6,574 | Common: 4,477
- Consistently up-regulated in both populations: 2,464
- Top candidates: SLC6A14, EPS8L3, PITX1, FXYD3, TNNT1, PVRL4 (NECTIN4), CEP55, NEK2

Plots and tables are in the `Result` folder (volcano plots, Venn diagram,
PCA, heatmap of top 30 genes, candidate gene list).

## Files
- `CCA_DE_analysis.R`: full analysis script
- `Result/`: output plots and tables

## Limitations
The two datasets come from different studies and platforms, so differences
between them are not purely population effects. Sample sizes are small.
Expression alone does not prove a gene is a vaccine target; these genes
are candidates for further validation.
