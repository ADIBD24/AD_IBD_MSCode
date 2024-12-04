# Load necessary libraries
library(edgeR)
library(affy)
library(oligo)
library(sva)
library(limma)
library(EnhancedVolcano)
library(TOmicsVis)
library(ggvenn)

# Step 1: Load RNA-seq datasets (example for ROSMAP, MSBB, and MayoRNAseq)
# Assume `rna_seq_counts` is a list of count matrices for RNA-seq datasets
tpm_transform <- function(count_matrix) {
  dge <- DGEList(counts = count_matrix)
  dge <- calcNormFactors(dge, method = "TMM")
  tpm <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
  return(tpm)
}

rna_seq_tpm <- lapply(rna_seq_counts, tpm_transform)

# Step 2: Load microarray datasets (example for GSE66407, GSE16879, and GSE179285)
# Assume `microarray_files` contains paths to CEL files for each dataset
rma_normalize <- function(file_paths) {
  raw_data <- read.celfiles(file_paths)
  eset <- rma(raw_data)
  return(exprs(eset))
}

microarray_normalized <- lapply(microarray_files, rma_normalize)

# Step 3: Combine RNA-seq and microarray data for batch effect correction
# Log2-transform RNA-seq TPM values
rna_seq_log2_tpm <- lapply(rna_seq_tpm, function(x) log2(x + 1))

# Combine all datasets
combined_data <- c(rna_seq_log2_tpm, microarray_normalized)

# Batch correction using ComBat
batch_labels <- c(rep("RNA-seq", length(rna_seq_log2_tpm)), rep("Microarray", length(microarray_normalized)))
combined_corrected <- ComBat(dat = do.call(cbind, combined_data), batch = batch_labels)

# Step 4: Perform PCA and validate batch correction
# PCA before correction
pca_before <- prcomp(do.call(cbind, combined_data))
plot(pca_before$x[, 1:2], col = as.factor(batch_labels), main = "PCA Before Correction")

# PCA after correction
pca_after <- prcomp(combined_corrected)
plot(pca_after$x[, 1:2], col = as.factor(batch_labels), main = "PCA After Correction")

# Step 5: Differential expression analysis with LIMMA
design <- model.matrix(~ 0 + factor(c("AD", "Control", "IBD")))
colnames(design) <- c("AD", "Control", "IBD")
fit <- lmFit(combined_corrected, design)
contrast.matrix <- makeContrasts(AD-Control, IBD-Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract DEGs with specified thresholds
deg_results <- topTable(fit2, coef = "AD-Control", number = Inf, adjust.method = "BH", p.value = 0.05, lfc = 0.2)

# Step 6: Visualize DEGs
# Volcano plot
EnhancedVolcano(deg_results,
                lab = rownames(deg_results),
                x = "logFC",
                y = "adj.P.Val",
                title = "Volcano Plot",
                pCutoff = 0.05,
                FCcutoff = 0.2)

# Step 7: Identify overlapping DEGs
# Create a list of DEGs from different datasets
deg_lists <- list(
  AD = rownames(deg_results[deg_results$logFC > 0.2 & deg_results$adj.P.Val < 0.05, ]),
  IBD = rownames(deg_results[deg_results$logFC < -0.2 & deg_results$adj.P.Val < 0.05, ])
)

# Venn diagram
ggvenn(deg_lists, fill_color = c("#FF9999", "#9999FF"))

# Step 8: Pathway enrichment and clustering (optional)
# Use TOmicsVis for cluster trends (example code)
cluster_trends <- clusterTrendPlot(deg_lists$AD, deg_lists$IBD)
print(cluster_trends)
