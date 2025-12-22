## MHT 

This repository supports the manuscript titled 'XXXXX' by providing the underlying methods and data analysis scripts to ensure reproducibility. As the manuscript is currently under peer review, the dataset provided here is simulated for demonstration purposes. This repository is intended specifically for reviewer evaluation. Upon formal publication, the complete and official dataset will be released.

## Details

### whole_clinic.wilcox.R: 
Longitudinal changes in the mKI and clinical indicators—including hormonal profiles, glucose metabolism, lipid profiles, inflammatory markers, coagulation parameters, complete blood count, and ultrasound measurements—were analyzed between baseline and post-intervention timepoints using the Wilcoxon signed-rank test for paired samples. mKI was evaluated at T4, T12, and T24, while clinical indicators were analyzed at T24 only. 

### whole.microbial.wilcox.R:
Differentially abundant taxa at each body site (vagina, urine, gut, tongue) during Menopausal Hormone Therapy (MHT), identified using the Wilcoxon signed-rank test between post-baseline timepoints (T04, T12, T24) to Baseline (BL).

### PREMANOVA.PHEN_TIME.R:
Permutational Multivariate Analysis of Variance (PERMANOVA) was used to assess whether overall microbial community structure (based on Bray-Curtis dissimilarity) was associated with phenotypic variation (e.g., mKI, sex hormones, blood glucose, blood lipids). Analyses were performed using the adonis2 function (vegan v2.7-1 in R).
		
### Spearman correlations.R: 
Pairwise Spearman correlations between delta values of phenotypic measures (e.g., ΔmKI, Δsex hormones, Δblood glucose, Δblood lipids, and other delta of continuous clinical indicators) were calculated using the cor.test function (stats v4.4.1 in R).
