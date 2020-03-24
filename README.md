# drug.gene.correlations

### Motivation
With the recent availability of drug sensitivity and gene vulnerability data across a wide representation of cancer cell lines, many attempts have been made to cross-analyze these datasets to elucidate cancer vulnerabilities and derive novel therapeutic targets. The most recent work in this area by [Gonçalves et al.](https://www.biorxiv.org/content/10.1101/2020.01.14.905729v1.full) applied linear mixed models to single drug-gene pairs to identify gene knockouts that most strongly correlated with the effects of specific drugs. Drug-gene pairs with statistically significant correlations were evaluated based on concordance with the drug’s putative target or proximity to it in protein interaction networks. Enrichment for on-target effects were found and novel biological insights were presented providing evidence for this method as a promising approach.

This data visualization is a continuation of Gonçalves et al.’s work extending their linear mixed model approach from single drug-gene pairs to regularized multiple regression of drugs with the entire gene set. By incorporating multiple regression to account for collinearities between gene data and LASSO regularization for feature selection, it is hypothesized that the resulting regression coefficients are better able to discriminate direct from indirect drug-gene relations. The ability for this model to account for multiple drug-gene interaction paradigms highlights its biological relevance and utility in screening new drugs and small-molecule compounds for therapeutic potential.

### Methods
Regularized multiple linear regression was done in R using the glmnet package with LASSO parameter λ chosen to minimize the mean-squared error during cross-validation. Single regression was done using drug-gene pair correlations.

Drug sensitivity data was provided by the Genomics of Drug Sensitivity in Cancer (GDSC) Project where 809 drugs were screened across 170 pan-cancer cell lines. Drug sensitivity was measured as the area-under-curve of the dose-response curve.

Gene vulnerability data was provided by Project Achilles from the Broad Institute where 18,334 genes were screened across 683 pan-cancer cell lines. Loss-of-function screen was done through CRISPR gene knockout. Gene vulnerability was measured as ratio of knockout cell viability to untreated control.

### Future Directions
While the dataset presented here can be used to highlight examples of where multiple regression outperforms single regression, a systematic analysis of the results still needs to be done. Future steps would entail construction of quantitative metrics to evaluate the validity of the model as well as functional testing to confirm whether results not corroborated by current knowledge are novel findings or simply artifacts due to unaccounted confounders or issues with the model. Once this is done, the validated model can then be leveraged to great effect in screening existing drugs for unknown effects and new compounds for therapeutic potential.
