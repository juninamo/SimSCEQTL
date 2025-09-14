
# SimSCEQTL: A Flexible Simulator for Single-Cell Genomics

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`SimSCEQTL` is a powerful and flexible R package for generating realistic single-cell RNA-seq datasets. Its primary purpose is to create complex, multi-individual datasets with a **known ground truth**, enabling rigorous benchmarking of computational methods for single-cell analysis.

The framework allows users to precisely control sources of variation, including cell type structure, batch effects, individual-level covariates, and cell-type-specific genetic effects (eQTLs).

---

## Key Features

`SimSCEQTL` provides a suite of functions to simulate each component of a single-cell experiment:

-   **🧬 Simulate Complex Cohorts:** Generate multi-individual metadata with customizable covariates such as age, sex, disease status, and batch.
-   **🔬 Model Cell Abundance:** Simulate realistic cell type compositions, including the ability to introduce differential abundance between groups (e.g., case vs. control) or complex interaction effects (e.g., disease effect is specific to a certain sex).
-   **📊 Generate Realistic Counts:** Use a latent feature model based on Gaussian Mixtures to generate sparse, overdispersed count data that mimics the properties of real scRNA-seq.
-   **🔗 Embed Ground-Truth eQTLs:** Simulate paired genotype data and embed cell-type-specific eQTL signals with user-defined effect sizes (correlation) and allele frequencies.
-   **⚙️ Enable Benchmarking of Analysis Pipelines:** The generated data is ideal for evaluating the statistical power and robustness of various analysis methods, such as pseudobulk approaches or cell-level mixed models.

---

## Installation

You can install the development version of `SimSCEQTL` from GitHub with:

```r
if (!require("remotes")) {
  install.packages("remotes")
}
remotes::install_github("juninamo/SimSCEQTL")
```

-----

## Example Workflow: Simulating a sc-eQTL Study

Here is a simplified example of how `SimSCEQTL` can be used to generate a dataset for benchmarking eQTL analysis pipelines.

```r
library(SimSCEQTL)
library(dplyr)
library(ggplot2)

# 1. Define simulation parameters
# Simulate 50 individuals, with cell type structure explaining 65% of variance
params <- list(
  n_individuals = 50,
  n_cells = 100,
  cluster_ratio = 0.65,
  batch_ratio = 0.1,
  n_thread = 4
)

# 2. Generate the complete dataset
# This single function call simulates metadata, latent features, counts, and genotypes.
sim_data <- generate_data(
  n_individuals = params$n_individuals,
  n_cells = params$n_cells,
  cluster_ratio = params$cluster_ratio,
  batch_ratio = params$batch_ratio,
  n_thread = params$n_thread,
  # Define a ground-truth eQTL in Cell Type 'A' for Gene1
  eqtl_celltype = "A",
  eqtl_gene = "Gene1",
  cor_vals = c(0, 0.3), # Simulate null (r=0) and active (r=0.3) eQTLs
  allele_freqs = c(0.1, 0.3)
)

# 3. Analyze the results (e.g., check power for pseudobulk eQTLs)
# The output contains all necessary data to perform power analysis,
# as demonstrated in our vignettes.
# For example, plot the ground-truth eQTL effect in the pseudobulk data.
genotype_df <- as.data.frame(t(sim_data$genotype_df_list$`0.3_0.3`))
colnames(genotype_df) <- paste0("SNP", 1:ncol(genotype_df))
genotype_df$subject_id <- rownames(genotype_df)

pseudobulk_data <- sim_data$bulk_df %>%
  filter(cell_type == "A") %>%
  left_join(genotype_df, by = "subject_id")

ggplot(pseudobulk_data, aes(x = as.factor(SNP1), y = Gene1)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(
    title = "Simulated eQTL in Cell Type A",
    subtitle = "Effect Size (r) = 0.3, MAF = 0.3",
    x = "Genotype Dosage",
    y = "Normalized Expression"
  ) +
  theme_minimal()
```

-----

## Vignettes

For detailed tutorials and use cases, please see our [vignettes](https://github.com/juninamo/SimSCEQTL/tree/master/vignettes):

  * **Simulating Cell Type Abundance:** A guide to simulating datasets with differential cell type proportions.
  * **Simulating Interaction Effects on Abundance:** Learn how to model complex interactions between covariates.
  * **Benchmarking sc-eQTL Methods:** A complete walkthrough of generating a dataset with ground-truth eQTLs and performing a power analysis.

-----

## 📝 Citation

If you use `SimSCEQTL` in your research, please cite our paper:

Jun Inamo, *et al*. A flexible R package for simulating single-cell RNA-seq data with ground-truth eQTLs for benchmarking analysis methods. (*XXXX*).

A preprint is available on *bioRxiv*: [https://doi.org/XXXX](https://www.google.com/search?q=https://doi.org/XXXX)

## Contact

For questions, bug reports, or issues, please contact:

**Name:** Jun Inamo  
**Email:** juninamo@keio.jp

Or open an issue on the [GitHub repository](https://www.google.com/search?q=https://github.com/juninamo/SimSCEQTL/issues).

  

## License

This repository is provided under the MIT License.



