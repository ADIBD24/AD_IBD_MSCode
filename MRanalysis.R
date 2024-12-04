# Load required packages
library(TwoSampleMR)
library(MRInstruments)

# Step 1: Load eQTL (exposure) and GWAS (outcome) data
# `eqtl_data` contains eQTL summary statistics for brain and intestinal tissues
# `gwas_data_AD` and `gwas_data_IBD` contain GWAS summary statistics for AD and IBD
eqtl_data <- read.table("eqtl_data.txt", header = TRUE)
gwas_data_AD <- read.table("gwas_data_AD.txt", header = TRUE)
gwas_data_IBD <- read.table("gwas_data_IBD.txt", header = TRUE)

# Step 2: Select instrumental variables (IVs)
# (1) Filter SNPs strongly associated with gene expression (P < 5×10^(-8))
eqtl_iv <- subset(eqtl_data, P < 5e-8)

# (2) Exclude SNPs associated with confounding traits using GWAS Catalog data
gwas_catalog <- read.table("gwas_catalog.txt", header = TRUE)
eqtl_iv <- eqtl_iv[!eqtl_iv$SNP %in% gwas_catalog$SNP, ]

# (3) Perform LD clumping to ensure independence among SNPs (r2 < 0.001, 10-Mb window)
eqtl_iv <- clump_data(eqtl_iv, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1e-8)

# Step 3: Harmonize data between exposure and outcome datasets
# Harmonize exposure (eQTL) with AD GWAS outcomes
harmonized_data_AD <- harmonise_data(exposure_dat = eqtl_iv, outcome_dat = gwas_data_AD)

# Harmonize exposure (eQTL) with IBD GWAS outcomes
harmonized_data_IBD <- harmonise_data(exposure_dat = eqtl_iv, outcome_dat = gwas_data_IBD)

# Step 4: Estimate causal effects using Wald ratio method
# For AD outcomes
mr_results_AD <- mr(harmonized_data_AD, method_list = c("mr_wald_ratio"))
mr_results_AD <- mr_results_AD[mr_results_AD$pval < 0.05, ] # Filter significant results

# For IBD outcomes
mr_results_IBD <- mr(harmonized_data_IBD, method_list = c("mr_wald_ratio"))
mr_results_IBD <- mr_results_IBD[mr_results_IBD$pval < 0.05, ] # Filter significant results

# Step 5: Extract and prioritize candidate genes
# Combine results from AD and IBD analyses
combined_results <- rbind(mr_results_AD, mr_results_IBD)

# Save the prioritized candidate genes for further analysis
write.table(combined_results, "prioritized_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Step 6: Validate IV assumptions (optional, diagnostic tests cannot be performed due to single IV per gene)
# This step is skipped due to limitations described in the study.

# Optional: Visualize results
# Scatter plot for AD outcomes
plot_mr_scatter(mr_results_AD, harmonized_data_AD)

# Scatter plot for IBD outcomes
plot_mr_scatter(mr_results_IBD, harmonized_data_IBD)

# Step 7: Interpret results
cat("Identified candidate genes with potential causal roles in AD and IBD are saved in 'prioritized_genes.txt'.\n")
