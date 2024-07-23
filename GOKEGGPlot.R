setwd("~/AD/IBD")

# Select differentially expressed genes (DEGs) that are either upregulated ("UP") or downregulated ("DOWN")
gene_entrezid <- DEG[DEG$g %in% c("UP", "DOWN"), ]
# Alternatively, use all DEGs for the analysis
gene_entrezid <- DEG

# Create a named vector of log fold changes (logFC) using entrez gene IDs as names
genelist <- gene_entrezid$logFC
names(genelist) <- gene_entrezid$ENTREZID

# Sort the genelist based on the log fold changes in a decreasing order
genelist <- sort(genelist, decreasing = TRUE)

# Perform GSEA analysis with the genelist
gsea_res <- gseGO(
  gene = genelist,
  OrgDb = "org.Hs.eg.db", # Organism database
  keyType = "ENTREZID", # Type of keys used in the gene list
  ont = "ALL", # Perform analysis for all GO ontologies
  pvalueCutoff = 1.0, # Cutoff for p-value
  pAdjustMethod = "BH", # Benjamini-Hochberg method for adjusting p-values
  minGSSize = 5, # Minimum gene set size
  maxGSSize = 500 # Maximum gene set size
)

# Print the first three results of the GSEA analysis
head(gsea_res, 3)

# Convert ENTREZID to gene symbols for readability
gsea_res_symbol <- setReadable(gsea_res, "org.Hs.eg.db", keyType = "ENTREZID")

# Plot the category network
cnetplot(
  gsea_res_symbol,
  showCategory = 3,
  color.params = list(foldChange = genelist) # Use genelist for color scaling
)

# Perform pairwise term similarity analysis
ora_pt <- pairwise_termsim(gsea_res_symbol)

# Plot the pairwise term similarity using a scatter plot
ssplot(
  ora_pt,
  showCategory = 30,
  drfun = NULL # Dimension reduction function
)

# Subset GSEA results to gene sets with a size less than 200
gsea_sub <- gsea_res[gsea_res$setSize < 200, ]

# Create a ridge plot for the enrichment distribution
ridgeplot(gsea_res) + labs(x = "enrichment distribution")

# Create a plot for the top 30 gene sets
gseaplot2(gsea_res, 1:30)

# Load the GseaVis library for advanced visualizations
library(GseaVis)

# Create a dot plot for the top 10 gene sets
dotplotGsea(
  data = gsea_res,
  topn = 10,
  pajust = 0.5
)

# Save the dot plot as a PDF file
ggsave(
  paste0(title, "-dotplotGseaGO.pdf"),
  width = 12,
  height = 8
)

# Create a volcano plot for the GSEA results
volcanoGsea(
  data = gsea_res,
  pajust = 0.5,
  p.adjust.CUTOFF = 0.5,
  nudge.y = c(-0.8, 0.8) # Adjust y-axis nudging for points
)

# Save the volcano plot as a PDF file
ggsave(
  paste0(title, "-volcanoGsea.pdf"),
  width = 10,
  height = 10
)


# Select top 3 mitochondrial-related gene sets for further analysis
# terms <- gsea_res[grepl("mitochondrial", gsea_res@result$Description),]$ID[1:3]
# terms <- gsea_res[grepl("immun", gsea_res@result$Description),]$ID[1:3] # Example for ROSMAP

# Perform GSEA network analysis for the selected terms
gseaNb(
  object = gsea_res,
  geneSetID = terms,
  subPlot = 2,
  termWidth = 35,
  legend.position = c(0.8, 0.8),
  addPval = T,
  pvalX = 0.05,
  pvalY = 0.05
)

# Save the network plot as a PDF file
ggsave(
  paste0(title, "-gseaNb-MT.pdf"),
  width = 10,
  height = 10
)

# Perform GO enrichment analysis for different ontologies (Biological Process, Cellular Component, Molecular Function)
ego_ALL <- enrichGO(
  gene = gene,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "BH",
  minGSSize = 1,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_CC <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_BP <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_MF <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)


# 4. Save the results to the current working directory
ego_ALL <- as.data.frame(ego_ALL) # Convert the GO enrichment results to a data frame
ego_result_BP <- as.data.frame(ego_BP) # Conversion for biological process results
ego_result_CC <- as.data.frame(ego_CC) # Conversion for cellular component results
ego_result_MF <- as.data.frame(ego_MF) # Conversion for molecular function results

# Combine the results of different GO ontologies into one data frame
ego <- rbind(ego_result_BP, ego_result_CC, ego_result_MF)

# Uncomment the following lines to write the results to CSV files
# write.csv(ego_ALL, file = "ego_ALL.csv", row.names = TRUE)
# write.csv(ego_result_BP, file = "ego_result_BP.csv", row.names = TRUE)
# write.csv(ego_result_CC, file = "ego_result_CC.csv", row.names = TRUE)
# write.csv(ego_result_MF, file = "ego_result_MF.csv", row.names = TRUE)
# write.csv(ego, file = "ego.csv", row.names = TRUE)

# 5. When there are too many enriched pathways, it's convenient to select a subset for visualization
display_number <- c(15, 15, 15) # The numbers represent the count of pathways to display for BP, CC, and MF

# Subset the GO enrichment results for the specified number of pathways
ego_result_BP <- ego_result_BP[1:display_number[1], ]
ego_result_CC <- ego_result_CC[1:display_number[2], ]
ego_result_MF <- ego_result_MF[1:display_number[3], ]

# Combine the selected pathways into a new data frame for GO enrichment
go_enrich_df <- data.frame(
  ID = c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
  Description = c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
  GeneNumber = c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type = factor(c(rep("biological process", display_number[1]),
                  rep("cellular component", display_number[2]),
                  rep("molecular function", display_number[3])),
                levels = c("biological process", "cellular component", "molecular function"))
)

# Shorten the pathway names to the first five words for better visualization
for (i in 1:nrow(go_enrich_df)) {
  description_splite <- strsplit(go_enrich_df$Description[i], split = " ")
  description_collapse <- paste(description_splite[[1]][1:5], collapse = " ")
  go_enrich_df$Description[i] <- description_collapse
  go_enrich_df$Description <- gsub(pattern = "NA", "", go_enrich_df$Description)
}

# Start plotting the GO bar charts
# Vertical bar chart
go_enrich_df$type_order <- factor(rev(as.integer(rownames(go_enrich_df))), labels = rev(go_enrich_df$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62") # Define colors for the bar chart

# Save the plot as a PDF file
pdf(file = paste0(paste0(title, gtitle), "-GO.pdf"))
ggplot(data = go_enrich_df, aes(x = type_order, y = GeneNumber, fill = type)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = COLS) +
  coord_flip() +
  xlab("GO term") +
  ylab("Gene_Number") +
  labs(title = "The Most Enriched GO Terms") +
  theme_bw()
dev.off()

# Horizontal bar chart
pdf(file = paste0(paste0(title, gtitle), '-GO.pdf'), width = 15, height = 8)
go_enrich_df$type_order <- factor(rev(as.integer(rownames(go_enrich_df))), labels = rev(go_enrich_df$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

ggplot(data = go_enrich_df, aes(x = type_order, y = GeneNumber, fill = type)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = COLS) +
  theme_bw() +
  xlab("GO term") +
  ylab("Num of Genes") +
  labs(title = "The Most Enriched GO Terms") +
  theme(axis.text.x = element_text(face = "bold", color = "gray50", angle = 80, vjust = 1, hjust = 1))
dev.off()

# Save the combined GO enrichment results as a CSV file
write.csv(ego, paste0(paste0(title, gtitle), '-GO.csv'))

# KEGG Enrichment
# Perform KEGG pathway enrichment analysis on a gene list
kk <- enrichKEGG(gene = gene, 
                 keyType = "kegg", 
                 organism = "human", 
                 qvalueCutoff = 0.05, 
                 pvalueCutoff = 0.05)

# Visualization
# Convert the enrichment results to a data frame for visualization
hh <- as.data.frame(kk)
# Set row names to integers and create an order for plotting
rownames(hh) <- 1:nrow(hh)
hh$order <- factor(rev(as.integer(rownames(hh))), labels = rev(hh$Description))

# Plot a bar chart to visualize the KEGG enrichment results
ggplot(hh, aes(y = order, x = Count, fill = p.adjust)) +
  geom_bar(stat = "identity", width = 0.7) +  # Set bar width
  #coord_flip() +  # Uncomment to flip the plot
  scale_fill_gradient(low = "red", high = "blue") +  # Set color gradient
  labs(title = "KEGG Pathways Enrichment",
       x = "Gene numbers", 
       y = "Pathways") +
  theme(axis.title.x = element_text(face = "bold", size = 16),
        axis.title.y = element_text(face = "bold", size = 16),
        legend.title = element_text(face = "bold", size = 16)) +
  theme_bw()
# Save the bar chart as a PDF file
ggsave(paste0(title, gtitle, "-enrichKK.pdf"))

# Plot a bubble chart to visualize the KEGG enrichment results
ggplot(hh, aes(y = order, x = Count)) +
  geom_point(aes(size = Count, color = -1 * p.adjust)) +  # Adjust point size and color
  scale_color_gradient(low = "green", high = "red") +
  labs(color = expression(p.adjust * size == "Count"),
       x = "Gene Number", 
       y = "Pathways",
       title = "KEGG Pathway Enrichment") +
  theme_bw()
# Save the bubble chart as a PDF file
ggsave(paste0(title, gtitle, "-enrichKK-1.pdf"))

# Write the KEGG enrichment results to a CSV file
write.csv(hh, paste0(title, gtitle, '-KEGG.csv'))

# Save the differentially expressed genes (DEG) data
save(gene_up, gene_down, DEG, file = paste0(title, "-DEG.Rdata"))
# Write the DEG data to a CSV file
write.csv(DEG, paste0(title, '-DEG.csv'))

# Load required libraries for plotting and analysis
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

# Define a function to annotate differentially expressed genes
deg_anno <- function(deg) {
  logFC_t <- 0.2  # Set threshold for log fold change
  # Annotate genes as UP, DOWN, or stable based on logFC and P.Value
  deg$g <- ifelse(deg$P.Value > 0.05, 'stable',
                  ifelse(deg$logFC > logFC_t, 'UP',
                         ifelse(deg$logFC < -logFC_t, 'DOWN', 'stable')))
  table(deg$g)  # Show the count of genes in each category
  head(deg)  # Show the first few rows of the DEG data
  deg$symbol <- rownames(deg)  # Set symbol as row names
  
  # Convert gene symbols to ENTREZID using org.Hs.eg.db
  df <- bitr(unique(deg$symbol), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  # Merge DEG data with the conversion result
  DEG <- merge(deg, df, by.y = 'SYMBOL', by.x = 'symbol')
  return(DEG)
}

# Apply the annotation function to DEG data from different conditions
TC_g <- deg_anno(TC_DEG)
HP_g <- deg_anno(HP_DEG)
EC_g <- deg_anno(EC_DEG)
FC_g <- deg_anno(FC_DEG)

# Filter out stable genes for Venn diagram analysis
TC_dg <- TC_g[TC_g$g != "stable",]
HP_dg <- HP_g[HP_g$g != "stable",]
EC_dg <- EC_g[EC_g$g != "stable",]
FC_dg <- FC_g[FC_g$g != "stable",]

# Draw a Venn diagram for the DEGs from different conditions
library("VennDiagram")
venn.diagram(x = list(TC_dg$symbol, HP_dg$symbol, EC_dg$symbol, FC_dg$symbol),
             scaled = F,  # Scaled by proportion
             alpha = 0.5,  # Transparency
             lwd = 1, lty = 1,  # Line width and type
             col = c('#9fc5e8', '#b6d7a8', "#ffe599", "#ead1dc"),  # Colors for each circle
             label.col = 'black',  # Color of the numbers inside the diagram
             cex = 2,  # Size of the numbers
             fontface = "bold",  # Font style
             fill = c('#9fc5e8', '#b6d7a8', "#ffe599", "#ead1dc"),  # Fill colors for the circles
             category.names = c("TC", "HP", "EC", "FC"),  # Names for each category
             cat.dist = c(0.2, 0.2, 0.1, 0.1),  # Distance of category labels from the diagram
             cat.pos = c(-20, 20, -20, 20),  # Angle of category labels
             cat.cex = 2,  # Size of category labels
             cat.fontface = "bold",  # Font style for category labels
             cat.col = c('#9fc5e8', '#b6d7a8', "#ffe599", "#ead1dc"),  # Color of category labels
             cat.default.pos = "outer",  # Default position for category labels
             output = TRUE,
             filename = './result/DEG/deg.png',  # File name for saving the diagram
             imagetype = "png",  # Image type
             resolution = 400,  # Resolution for the image
             compression = "lzw")  # Compression method

# Draw the Venn diagram using grid.draw
grid.draw(data)

