# Set the working directory to the specified path
setwd("~/AD/IBD/")

# Load necessary R libraries
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(OmnipathR)

# Load the processed RNA-seq data and annotations for different datasets
load("~/AD/IBD/result/DEG/GSE66407-deg.Rdata")
# ... (other load commands for different datasets)

# Define a title for the dataset being analyzed
title = "GSE66407"

# Convert the explan_final data frame to a data frame and remove NA values, setting row names
counts <- as.data.frame(explan_final)

# Create a design data frame from sample information and group metadata, and remove NA values
design <- data.frame(Sample = rownames(pdata), Group = pdata$Group)

# Extract t-values per gene from the deg dataset and set row names
degg <- data.frame(t=deg$t)
rownames(degg) <- rownames(deg)

# Retrieve the PROGENy network for human, which contains curated pathways and their target gene sets
net <- get_progeny(organism = 'human', top = 500)

# Run the multivariate linear model (mlm) to infer pathway activities from gene expression data
sample_acts <- run_mlm(mat=counts, net=net, .source='source', .target='target', .mor='weight', minsize = 5)

# Transform the sample activities to a wide matrix format for visualization
sample_acts_mat <- sample_acts %>%
  pivot_wider(id_cols = 'condition', names_from = 'source', values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

# Scale the sample activities matrix per feature and choose a color palette for the heatmap
sample_acts_mat <- scale(sample_acts_mat)
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white", "red"))(palette_length)

# Plot the heatmap of sample activities and export the PDF file
pdf(paste0(title,"-pheatmap-decoupleR.pdf"),width = 10,height = 12)
pheatmap(sample_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 
dev.off()

# Run the mlm again with the degg dataset to contrast pathway activities between conditions
contrast_acts <- run_mlm(mat=degg, net=net, .source='source', .target='target', .mor='weight', minsize = 5)

# Plot the contrast of pathway activities and export the PDF file
pdf(paste0(title,"-DEG-Pathways-decoupleR.pdf"),width = 10,height = 8)
ggplot(contrast_acts, aes(x = reorder(source, score), y = score)) + 
  # ... (plotting code with theme and axis settings)
  dev.off()

# For a specific pathway (e.g., 'JAK-STAT'), filter the network, extract relevant genes, and visualize their activities
pathway <- 'JAK-STAT'
df <- net %>%
  filter(source == pathway) %>%
  # ... (data manipulation and plotting code)
  write.csv(df,paste0(title,"-DEG-",pathway,"-decoupleR-gene.csv"),row.names = F)

# Plot the activities of the filtered genes and export the PDF file
pdf(paste0(title,"-DEG-",pathway,"-decoupleR.pdf"),width = 10,height = 8)
ggplot(df, aes(x = weight, y = t_value, color = color)) +
  # ... (plotting code with geom and theme settings)
  dev.off()

# Use the CollecTRI network to infer transcription factor activities from gene expression data
net <- get_collectri(organism='human', split_complexes=FALSE)

# Run the univariate linear model (ulm) to infer transcription factor activities
sample_acts <- run_ulm(mat=counts, net=net, .source='source', .target='target', .mor='mor', minsize = 5)

# Visualize the top transcription factors with the most variable activities across samples
n_tfs <- 25
# ... (data transformation and plotting code for the heatmap)

# Contrast the activities of transcription factors between conditions
contrast_acts <- run_ulm(mat=degg[, 't', drop=FALSE], net=net, .source='source', .target='target', .mor='mor', minsize = 5)
# ... (filtering and plotting code for the transcription factor activities)

# Visualize the target genes of a specific transcription factor (e.g., 'RELA') and their association with the factor's activity
tf <- 'RELA'
df <- net %>%
  filter(source == tf) %>%
  # ... (data manipulation and plotting code)
  write.csv(df,paste0(title,"-DEG-",tf,"-decoupleR-gene.csv"),row.names = F)

# Plot the genes associated with the transcription factor and export the PDF file
pdf(paste0(title,"-DEG-",tf,"-decoupleR.pdf"),width = 10,height = 8)
ggplot(df, aes(x = logfc, y = -log10(p_value), color = color, size=abs(mor))) +
  # ... (plotting code with geom and theme settings)
  dev.off()