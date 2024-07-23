library(WGCNA)
setwd("~/AD/IBD/wgcna")

### Draw a hierarchical clustering dendrogram for samples
if (T) {
  # Perform hierarchical clustering on samples based on expression data
  sampleTree <- hclust(dist(datExpr0), method = "average")
  # Plot the dendrogram
  par(mar = c(0, 5, 2, 0))  # Set the margins of the plot
  plot(sampleTree, main = "Sample clustering", sub = "", xlab = "", cex.lab = 2, 
       cex.axis = 1, cex.main = 1, cex.lab = 1)
  
  # If samples have traits or phenotypes, add corresponding colors to check clustering合理性
  sample_colors <- numbers2colors(as.numeric(factor(datTraits$Group)), 
                                  colors = rainbow(length(table(datTraits$Group))), 
                                  signed = FALSE)
  # Plot the dendrogram with sample traits
  par(mar = c(1, 4, 3, 1), cex = 0.8)
  pdf(paste0(title, "-step1_Sample dendrogram and trait-rm.pdf"), width = 10, height = 6)
  plotDendroAndColors(sampleTree, sample_colors,
                      groupLabels = "trait",
                      cex.dendroLabels = 0.8,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.01,
                      main = "Sample dendrogram and trait")
  # Add a line to show the cut
  abline(h = 110, col = "red")  # Adjust the height based on actual situation
  dev.off()
}

### Remove significant outliers if they exist
if (T) {
  # Cut the dendrogram to remove outliers based on height and minimum cluster size
  clust <- cutreeStatic(sampleTree, cutHeight = 150, minSize = 10)  # Adjust cutHeight based on actual situation
  table(clust)
  # Subset the data to keep only the desired samples
  keepSamples <- (clust == 1)
  datExpr0 <- datExpr0[keepSamples, ]
  datTraits <- datTraits[keepSamples, ]
  datTraits <- datTraits[complete.cases(datTraits), ]
  dim(datExpr0)  # Check the dimensions of the expression data
}

### Assess data quality with PCA
library(factoextra)
library(FactoMineR)
# Perform PCA on the expression data
dat.pca <- PCA(datExpr0, graph = F)
# Visualize the PCA results with sample grouping
pca <- fviz_pca_ind(dat.pca,
                    title = "Principal Component Analysis",
                    legend.title = "Groups",
                    geom.ind = c("point", "text"),  # "point", "text"
                    pointsize = 2,
                    labelsize = 4,
                    repel = TRUE,  # Avoid overlapping labels
                    col.ind = group_list,  # Color by group
                    axes.linetype = NA)  # Remove axes lines
pca
# Save the PCA plot as a PDF file
ggsave(pca, filename = "mayo-cqn-step1_Sample PCA analysis.pdf", width = 15, height = 15)

# Set the R-squared cut-off for network topology analysis
R.sq_cutoff = 0.8
if (T) {
  # Perform network topology analysis with different powers
  powers <- c(seq(1, 20, by = 1), seq(22, 30, by = 2))
  sft <- pickSoftThreshold(datExpr0, 
                           networkType = "unsigned",
                           powerVector = powers, 
                           RsquaredCut = R.sq_cutoff,  
                           verbose = 5)
  # Plot the results to find the optimal power value
  pdf(paste0(title, "-step2_power-value.pdf"), width = 16, height = 12)
  par(mfrow = c(1, 2))
  cex1 = 0.9
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n")
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, cex = cex1, col = "red")
  abline(h = R.sq_cutoff, col = "red")
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n")
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
  abline(h = 100, col = "red")
  dev.off()
}

# Determine the optimal power value for the network
power = sft$powerEstimate

# If no suitable power value is found, adjust based on sample count and network type
if (is.na(power)) {
  type = "unsigned"
  nSamples = nrow(datExpr)
  power = ifelse(nSamples < 20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples < 30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples < 40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))))
}

# Save the input data for further analysis
save(nGenes, nSamples, datExpr0, datTraits, design, file = paste0(title, "-step1_input.Rdata"))

###################### 3. Construct Weighted Co-expression Network and Identify Gene Modules ####################

# Load previously saved data
load(file = "step1_input.Rdata")
load(file = "step2_power_value.Rdata")

# Define the correlation function to be used
cor <- WGCNA::cor

if (T) {
  # Construct gene modules using blockwiseModules
  net <- blockwiseModules(
    datExpr0,  # Expression data
    power = power,  # Power value from previous step
    maxBlockSize = ncol(datExpr0),  # Maximum block size for computation
    corType = "pearson",  # Correlation type
    networkType = "unsigned",  # Network type
    TOMType = "unsigned",  # TOM matrix type
    minModuleSize = 30,  # Minimum module size
    mergeCutHeight = 0.25,  # Merge modules with higher threshold
    numericLabels = TRUE,  # Return numeric labels for modules
    saveTOMs = F,  # Save TOMs for future use
    verbose = 3,
    nThreads = 3  # Number of threads for parallel computation
  )
  # Show the count of genes in each module
  table(net$colors)
}

# Visualization of modules using hierarchical clustering dendrogram
if (T) {
  # Convert module labels to colors
  moduleColors <- labels2colors(net$colors)
  table(moduleColors)
  
  # Plot the dendrogram and module colors
  pdf(paste0(title, "-step3_genes-modules_ClusterDendrogram.pdf"), width = 16, height = 12)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}

# Save the network and module colors for future use
save(net, moduleColors, file = paste0(title, "-step3_genes_modules.Rdata"))

##################### Distribution-based Network Construction, generally not used #################################
if (T) {
  # Construct a weighted co-expression network in two steps:
  # 1. Calculate the adjacency matrix (correlation coefficients between gene expressions)
  # 2. Convert the adjacency matrix to a topological overlap matrix (TOM)
  
  # (1) Network construction: Co-expression similarity and adjacency
  adjacency = adjacency(datExpr0, power = power)
  
  # (2) Convert adjacency to topological overlap
  TOM = TOMsimilarity(adjacency)
  dissTOM = 1 - TOM
  
  # (3) Cluster the TOM using hierarchical clustering
  geneTree = hclust(as.dist(dissTOM), method = "average")
  # Plot the gene clustering dendrogram
  plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04)
  
  # (4) Dynamic tree cut for module identification
  # Set the minimum module size
  minModuleSize = 20
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  table(dynamicMods)
  
  # Convert numeric labels to colors
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  
  # Plot the gene dendrogram and module colors
  pdf(file = paste0(title, "-Gene dendrogram and module colors.pdf"), width = 16, height = 12)
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  dev.off()
  
  # (5) Merge modules with similar expression profiles
  # Calculate eigengenes for each module
  MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1 - cor(MEs)
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average")
  # Set the merging height cut (e.g., 0.25 for 75% similarity)
  MEDissThres = 0.25
  # Plot the module eigengene clustering
  pdf(file = paste0(title, "-Clustering of module eigengenes.pdf"), width = 16, height = 12)
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  abline(h = MEDissThres, col = "red")
  dev.off()
  
  # Merge close modules automatically
  merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors
  # Eigengenes of the new merged modules
  mergedMEs = merge$newMEs
  # Count of genes in each merged module
  table(mergedColors)
  
  # Plot the gene dendrogram with dynamic and merged module colors
  pdf(file = paste0(title, "-step3_stepbystep_DynamicTreeCut_genes-modules.pdf"), width = 16, height = 12)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
  # Save the step-by-step gene modules data
  moduleColors = mergedColors
  colorOrder = c("grey", standardColors(50))
  moduleLabels = match(moduleColors, colorOrder) - 1
  MEs = mergedMEs
  save(MEs, mergedColors, moduleLabels, moduleColors, geneTree,
       file = paste0(title, "-step3_stepByStep_genes_modules.Rdata"))
}

####################### 4. Associate gene modules with phenotypes #####################################

# Clear the current R environment
rm(list = ls())

# Load previously saved data
load(file = "step1_input.Rdata")
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")

# Save the number of genes and samples, along with the expression data and traits
nGenes <- ncol(datExpr0)
nSamples <- nrow(datExpr0)
save(nGenes, nSamples, datExpr0, datTraits, design, file = paste0(title, "-step1_input.Rdata"))

## Module-trait correlation heatmap
if (T) {
  # Convert traits to factors for the analysis
  datTraits$Group <- as.factor(datTraits$Group)
  datTraits$`tissue:ch1` <- as.factor(datTraits$`tissue:ch1`)
  datTraits$`inflammation:ch1` <- as.factor(datTraits$`inflammation:ch1`)
  
  # Create the design matrix for the analysis
  design <- model.matrix(~0 + datTraits$Group + datTraits$`tissue:ch1` + datTraits$`inflammation:ch1`)
  
  # Calculate module eigengenes and order them
  MES0 <- moduleEigengenes(datExpr0, moduleColors)$eigengenes
  MEs <- orderMEs(MES0)
  
  # Compute the correlation between module eigengenes and traits
  moduleTraitCor <- cor(MEs, design, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  
  # Create a text matrix for the correlation values and p-values
  textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  # Plot the module-trait correlation heatmap
  pdf(file = paste0(title, "-step4_Module-trait-relationship_heatmap.pdf"),
      width = 2 * length(colnames(design)),
      height = 0.6 * length(names(MEs)))
  par(mar = c(5, 9, 3, 3))  # Set the margins for the plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-0.6, 0.6),
                 main = "Module-trait relationships")
  dev.off()
  save(design, file = paste0(title, "-step4_design.Rdata"))
}

## Module-trait correlation boxplot
if (T) {
  # Merge module eigengenes with trait data for further analysis
  mes_group <- merge(MEs, datTraits, by = "row.names")
  
  # Draw boxplots for the module-trait relationship
  draw_ggboxplot <- function(data, Module = "Module", group = "group") {
    ggboxplot(data, x = group, y = Module,
              ylab = paste0(Module),
              xlab = group,
              fill = group,
              palette = "jco",
              legend = "") + stat_compare_means()
  }
  
  # Batch plot boxplots for each module
  colorNames <- names(MEs)
  pdf(paste0(title, "-step4_Module-trait-relationship_boxplot-Disease.pdf"), width = 25, height = 1.6 * ncol(MEs))
  p <- lapply(colorNames, function(x) {
    draw_ggboxplot(mes_group, Module = x, group = "Disease")
  })
  do.call(grid.arrange, c(p, ncol = 4))  # Arrange plots in a grid
  dev.off()
}

## Gene-module and gene-trait significance scatter plot
# Genes that are significantly associated with both a module and a trait may be particularly important.

# Select discrete traits for analysis
choose_group <- "AD"

if (T) {
  # Calculate the module membership correlation matrix
  geneModuleMembership <- as.data.frame(cor(datExpr0, MEs, use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) <- paste0("MM", modNames)
  names(MMPvalue) <- paste0("p.MM", modNames)
  
  # Calculate the gene significance matrix for traits
  trait <- as.data.frame(design[, choose_group])
  geneTraitSignificance <- as.data.frame(cor(datExpr0, trait, use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) <- paste0("GS")
  names(GSPvalue) <- paste0("GS")
  
  # Visualize gene-module and gene-trait significance
  selectModule <- modNames  # Select all modules for plotting
  pdf(paste0(title, "-step4_gene-Module-trait-significance.pdf"), width = 15, height = 1.5 * ncol(MEs))
  par(mfrow = c(ceiling(length(selectModule) / 3), 3))  # Set up a matrix for multiple plots
  for (module in selectModule) {
    column <- match(module, selectModule)
    print(module)
    moduleGenes <- moduleColors == module
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = "Gene significance for trait",
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  }
  dev.off()
}

########################## 5. WGCNA Visualization: TOMplot and Eigengene-adjacency-heatmap ##################################

if (T) {
  # Compute the TOM (Topological Overlap Matrix) from expression data
  TOM = TOMsimilarityFromExpr(datExpr0, power = power)
  # Calculate the dissimilarity matrix
  dissTOM = 1 - TOM
  
  # Draw the TOMplot for all genes
  if (T) {
    # Obtain the gene dendrogram
    geneTree = net$dendrograms[[1]]
    # Raise the power to enhance contrast
    plotTOM = dissTOM^7
    # Set the diagonal to NA
    diag(plotTOM) = NA
    # Save the TOMplot as a PNG image
    png(paste0(title, "-step5_TOMplot_Network-heatmap.png"), width = 800, height = 600)
    # Draw the TOMplot
    TOMplot(plotTOM, geneTree, moduleColors,
            col = gplots::colorpanel(250, 'red', "orange", 'lemonchiffon'),
            main = "Network heatmap plot")
    dev.off()
  }
  
  # Draw the TOMplot for selected genes to save time (optional)
  if (F) {
    # Select a subset of genes
    nSelect = 0.1 * nGenes
    set.seed(123)
    select = sample(nGenes, size = nSelect)
    # Compute the dissimilarity matrix for selected genes
    selectTOM = dissTOM[select, select]
    # Obtain the dendrogram for selected genes
    selectTree = hclust(as.dist(selectTOM), method = "average")
    # Extract module colors for selected genes
    selectColors = moduleColors[select]
    # Raise the power to enhance contrast
    plotDiss = selectTOM^7
    # Set the diagonal to NA
    diag(plotDiss) = NA
    # Save the TOMplot as a PDF document
    pdf(paste0(title, "step5_select_TOMplot_Network-heatmap.pdf"), width = 8, height = 6)
    # Draw the TOMplot for selected genes
    TOMplot(plotDiss, selectTree, selectColors,
            col = gplots::colorpanel(250, 'red', "orange", 'lemonchiffon'),
            main = "Network heatmap plot of selected genes")
    dev.off()
  }
}

### Module correlation display: Eigengene-adjacency-heatmap
if (T) {
  # Calculate module eigengenes
  MEs = moduleEigengenes(datExpr0, moduleColors)$eigengenes
  # Order the eigengenes
  MET = orderMEs(MEs)
  
  # If adding phenotype data
  if (T) {
    # Continuous trait
    #MET = orderMEs(cbind(MEs, datTraits$Age))
    # Categorical trait, convert to 0-1 matrix (stored in design)
    #AD = as.data.frame(design[,1])
    #names(AD) = "AD"
    # Add the weight to existing module eigengenes
    #MET = orderMEs(cbind(MEs, AD))
  }
  # Save the Eigengene-adjacency-heatmap as a PDF document
  pdf(paste0(title, "-step5_module_cor_Eigengene-dendrogram-AD.pdf"), width = 10, height = 12)
  # Plot the eigengene network
  plotEigengeneNetworks(MET, setLabels = "",
                        marDendro = c(0, 4, 1, 4),  # Margins for the dendrogram
                        marHeatmap = c(5, 5, 1, 2), # Margins for the heatmap
                        cex.lab = 0.8,
                        xLabelsAngle = 90)
  dev.off()
}

#################### 6. Select gene modules of interest for GO analysis ####################
### Condition settings
OrgDb = "org.Hs.eg.db"  # "org.Mm.eg.db"  "org.Hs.eg.db"
genetype = "SYMBOL"    # "SYMBOL"   "ENSEMBL"

# List available modules
table(moduleColors)

# Select modules of interest
choose_module <- c("violet", "royalblue", "green", "darkmagenta", "darkorange2", "")
for (i in choose_module) {
  module = i
  if (T) {
    # Load required libraries for GO analysis
    library(clusterProfiler)
    library(org.Hs.eg.db)
    
    # Prepare the gene module data frame
    gene_module <- data.frame(gene = colnames(datExpr0),
                              module = moduleColors)
    # Save the gene module data frame as a CSV file
    write.csv(gene_module, paste0(title, file = "step6_gene_moduleColors.csv"), row.names = F, quote = F) 
    
    # Convert gene identifiers to ENTREZID
    tmp <- bitr(gene_module$gene, fromType = genetype,
                toType = "ENTREZID",
                OrgDb = OrgDb)
    # Merge the data frames
    gene_module_entrz <- merge(tmp, gene_module, by.x = genetype, by.y = "gene")
    
    # Select genes belonging to the module of interest
    choose_gene_module_entrz <- gene_module_entrz[gene_module_entrz$module %in% module,]
    
    ### Run GO analysis
    formula_res <- compareCluster(
      ENTREZID ~ module,
      data = choose_gene_module_entrz,
      fun = "enrichGO",
      OrgDb = OrgDb,
      ont = "BP",  # One of "BP", "MF", and "CC" or "ALL"
      pAdjustMethod = "BH",
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.1
    )
    # Save the GO analysis results as a CSV file
    write.csv(formula_res@compareClusterResult,
              file = paste0(paste(title, module, sep = "-"), "-step6_module_GO_term.csv"))
    
    # Simplify the GO analysis results
    lineage1_ego <- clusterProfiler::simplify(
      formula_res,
      cutoff = 0.25,
      by = "p.adjust",
      select_fun = min
    )
    # Save the simplified GO analysis results
    save(gene_module, formula_res, lineage1_ego, file = paste0(paste(title, module, sep = "-"), "step6_module_GO_term.Rdata"))
    write.csv(lineage1_ego@compareClusterResult,
              file = paste0(paste(title, module, sep = "-"), "-step6_module_GO_term2.csv"))
    
    ### Plot a dotplot for GO analysis results
    dotp <- dotplot(formula_res,
                    showCategory = 10,
                    includeAll = TRUE, # Show results with overlap
                    label_format = 90)
    ggsave(dotp, filename = paste0(paste(title, module, sep = "-"), "-step6_module_GO_term.pdf"),
           width = 15,
           height = 15)
  }
}

############################### 7. Visualize genes in the gene module of interest as a heatmap ######################################

# Clear the current R environment
rm(list = ls())

# Load previously saved data
load(file = 'step1_input.Rdata')
load(file = "step3_genes_modules.Rdata")

# List available modules and their colors
table(moduleColors)

# Set the title and module of interest
title = "TC"
module = "green"

### Visualize the gene module of interest as a heatmap
if (T) {
  # Subset the data for the module of interest
  dat = datExpr0[, moduleColors == module]
  
  # Load the pheatmap library for heatmap visualization
  library(pheatmap)
  
  # Scale the data and transpose the matrix for visualization
  n = t(scale(dat))
  
  # Define group list for annotation
  group_list = datTraits$Disease
  ac = data.frame(g = group_list)
  rownames(ac) = colnames(n)
  
  # Generate the heatmap with annotations
  pheatmap::pheatmap(n,
                     fontsize = 8,
                     show_colnames = T,
                     show_rownames = F,
                     cluster_cols = T,
                     annotation_col = ac,
                     width = 8,
                     height = 6,
                     angle_col = 45,
                     main = paste0("module_", module, "-gene heatmap"),
                     filename = paste0(title, "-step7_module_", module, "_Gene-heatmap.pdf"))
}

################### 8. Export genes in the module of interest for network analysis in VisANT or Cytoscape ######################

# Compute the TOM (Topological Overlap Matrix) from expression data
TOM <- TOMsimilarityFromExpr(datExpr0, power = power)

# Save the TOM for future use
save(TOM,
     file = paste0(title, "-step8_TOM.Rdata"))

# Select a module of interest
choose_module = "violet"

# Loop through the selected modules
for (i in choose_module) {
  module = i
  if (T) {
    # Extract gene names for the module of interest
    gene <- colnames(datExpr0)
    inModule <- moduleColors == module
    modgene <- gene[inModule]
    
    # Subset the TOM for the module of interest
    modTOM <- TOM[inModule, inModule]
    dimnames(modTOM) <- list(modgene, modgene)
    
    # Select top connected genes to focus on the most relevant ones
    nTop = 100
    IMConn = softConnectivity(datExpr0[, modgene]) # Compute connectivity
    top = (rank(-IMConn) <= nTop) # Select top connected genes
    filter_modTOM <- modTOM[top, top]
    
    # Export the network for VisANT
    vis <- exportNetworkToVisANT(filter_modTOM,
                                 file = paste("step8_visANTinput-", module, ".txt", sep = ""),
                                 weighted = T, threshold = 0)
    
    # Export the network for Cytoscape
    cyt <- exportNetworkToCytoscape(filter_modTOM,
                                    edgeFile = paste(paste0(title, "-step8_CytoscapeInput-edges-"), paste(module, collapse = "-"), ".txt", sep = ""),
                                    nodeFile = paste(paste0(title, "-step8_CytoscapeInput-nodes-", paste(module, collapse = "-"), ".txt", sep = "")),
                                    weighted = TRUE,
                                    threshold = mean(vis$weight), # Weighted threshold for selection
                                    nodeNames = modgene[top],
                                    nodeAttr = moduleColors[inModule][top])
  }
}

# Add annotations for MT-related nodes in the network
nodeFile = read.table(paste(paste0(title, "-step8_CytoscapeInput-nodes-", paste(module, collapse = "-"), ".txt", sep = "")), header = 0, skip = 1)
# Define annotations based on gene categories
for (i in rownames(nodeFile)) {
  if (nodeFile[i, "V1"] %in% mt_gene) {
    nodeFile[i, "Type"] = "MT-located"
  } else if (nodeFile[i, "V1"] %in% ep_gene$V2) {
    nodeFile[i, "Type"] = "MT-epstasis"
  } else if (nodeFile[i, "V1"] %in% intersect(ep_gene$V2, mt_gene)) {
    nodeFile[i, "Type"] = "MT-located and epstasis"
  } else {
    nodeFile[i, "Type"] = "Other"
  }
}
# Save the annotated node file
write.csv(nodeFile, paste(paste0(title, "-step8_CytoscapeInput-nodes-"), paste(module, collapse = "-"), ".txt", sep = ""), row.names = FALSE)
##添加网络中MT相关节点的标注
nodeFile = read.table(paste(paste0(title,"-step8_CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep="")),header =0,skip = 1)
for(i in rownames(nodeFile)){
  if(nodeFile[i,]$V1 %in% mt_gene){
    nodeFile[i,"Type"]="MT-located"}
  else if (nodeFile[i,]$V1 %in% ep_gene$V2){
    nodeFile[i,"Type"]="MT-epstasis"}
  else if (nodeFile[i,]$V1 %in% intersect(ep_gene$V2,mt_gene)){
    nodeFile[i,"Type"]="MT-located and epstasis"}
  else{ nodeFile[i,"Type"]="Other"}
}
write.csv(nodeFile,paste(paste0(title,"-step8_CytoscapeInput-nodes-"), paste(module, collapse="-"), ".txt", sep=""),row.names = FALSE)


##########################
##线粒体相关基因富集
#https://www.jianshu.com/p/360dc479c685
mt_gene = read.csv("/home/xxu/AD/The RNAseq Harmonization Study/MT-gene.csv", header=0)
ep_gene = read.csv("/home/xxu/AD/The RNAseq Harmonization Study/intersnp-com-symbol.csv", header=0)
mt_gene <- unique(mt_gene[mt_gene$V7>=3,]$V2)

title = "GSE179285"
setwd("~/AD/IBD/result/wgcna/module/")

# Load the WGCNA data and module information
load("~/AD/IBD/result/wgcna/GSE179285-step1_input.Rdata")
load("~/AD/IBD/result/wgcna/GSE179285-step8_TOM.Rdata")
load("~/AD/IBD/result/wgcna/GSE179285-step3_stepByStep_genes_modules.Rdata")

# Loop through each unique module color to process the enrichment of mitochondrial-related genes
for (i in unique(moduleColors)){
  module = i
  cat("PROCESS:",title ,"\n")  
  if(T){
    ### Extract gene names of the module of interest
    gene <- colnames(datExpr0) 
    inModule <- moduleColors==module
    modgene <- gene[inModule]
    
    # Calculate the enrichment of mitochondrial genes in the module using hypergeometric test
    a = length(intersect(modgene,mt_gene))
    p1 <- phyper(a-1,length(mt_gene),length(gene)-length(mt_gene),length(modgene),lower.tail=F)
    cat(module,"-MT-located gene phyper:",p1,"\n")
    
    # Calculate the enrichment of epistatic genes in the module using hypergeometric test
    b = length(intersect(modgene,ep_gene$V2))
    p2 <- phyper(b-1,length(ep_gene$V2),length(gene)-length(ep_gene$V2),length(modgene),lower.tail=F)
    cat(module,"-MT-epstasis gene phyper:",p2,"\n")
    
    # If the p-value is significant, calculate the KME values and identify hub genes
    if (p1<=0.05 | p2<=0.05){
      datKME=signedKME(datExpr0, MEs, outputColumnName="kME_MM.")
      # ... (rest of the code for identifying hub genes)
    }
  }
}


#####Heatmap with the proportion of mitochondrial-related genes

datTraits$Group <- as.factor(datTraits$Group)
datTraits$`tissue:ch1` = as.factor(datTraits$`tissue:ch1`)
datTraits$`inflammation:ch1` = as.factor(datTraits$`inflammation:ch1`)
design <- model.matrix(~0+datTraits$Group+datTraits$`tissue:ch1`+datTraits$`inflammation:ch1`)
#  colnames(design) <- c("AD","Control","Pathologic Aging","PSP","Sex","Age At Death","APOE4-1","APOE4-2") #get the group
#colnames(design) <- c("AD","Control","Sex","Age") #get the group
colnames(design) <- c("Control","IBD","Tissue:Ileum","Inflammation") #get the group

MES0 <- moduleEigengenes(datExpr0,moduleColors)$eigengenes  #Calculate module eigengenes.
MEs <- orderMEs(MES0)  #Put close eigenvectors next to each other

# Calculate the correlation between module eigengenes and traits, and the p-values
moduleTraitCor <- cor(MEs,design,use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)

# Prepare the text matrix for the heatmap with the proportions and p-values of mitochondrial-related genes
textMatrix <- paste0(signif(moduleTraitCor,2),"\n(",
                     signif(moduleTraitPvalue,2),")")
dim(textMatrix) <- dim(moduleTraitCor)
textMatrix <- as.data.frame(textMatrix)
colnames(textMatrix) <- colnames(design)
rownames(textMatrix) <- colnames(MEs)
moduleTraitCor <- as.data.frame(moduleTraitCor)

moduleTraitCor$"MT-located" <- rep(1,length(colnames(MEs)))
moduleTraitCor$"MT-epistasis" <- rep(1,length(colnames(MEs)))
textMatrix$"MT-located" <- rep(1,length(colnames(MEs)))
textMatrix$"MT-epistasis" <- rep(1,length(colnames(MEs)))

for (i in unique(moduleColors)){
  module = i
  cat("PROCESS:",title ,"\n")  
  if(T){
    ### select interesting module
    gene <- colnames(datExpr0) 
    inModule <- moduleColors==module
    modgene <- gene[inModule]
    
    a = length(intersect(modgene,mt_gene))
    per1 <- a/length(modgene)
    p1 <- phyper(a-1,length(mt_gene),length(gene)-length(mt_gene),length(modgene),lower.tail=F)
    cat(module,"-MT-located gene phyper:",per1," ",p1,"\n")
    textMatrix[paste0("ME",module),"MT-located"] = paste0(signif(per1,2),"\n(",
                                                          signif(p1,2),")")
    moduleTraitCor[paste0("ME",module),"MT-located"]=per1
    
    b = length(intersect(modgene,ep_gene$V2))
    per2 <- b/length(modgene)
    p2 <- phyper(b-1,length(ep_gene$V2),length(gene)-length(ep_gene$V2),length(modgene),lower.tail=F)
    cat(module,"-MT-epstasis gene phyper:",per2," ",p2,"\n")
    textMatrix[paste0("ME",module),"MT-epistasis"] = paste0(signif(per2,2),"\n(",
                                                            signif(p2,2),")")
    moduleTraitCor[paste0("ME",module),"MT-epistasis"]=per2
  }  }

pdf(file = paste0(title,"-Module-trait-MT-relationship_heatmap.pdf"),
    width = 2*(length(colnames(design))+2), 
    height = 0.6*length(names(MEs)) )
par(mar=c(5, 9, 3, 3)) 
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(textMatrix),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = F,
               colors = colorRampPalette(brewer.pal(11, "BrBG"))(50),
               textMatrix = textMatrix,
               setStdMargins = F,
               cex.text = 0.5,
               zlim = c(-1.0,1.0), 
               main = paste0(title,"-Module-trait relationships"))
dev.off()

for (i in rownames(moduleTraitCor)){
  moduleTraitCor[i,"rank"] <- sum(abs(moduleTraitCor[i,"IBD"]),mean(moduleTraitCor[i,"MT-located"],moduleTraitCor[i,"MT-epistasis"]))
}
moduleTraitCor <- moduleTraitCor[order(moduleTraitCor$rank,decreasing = T),]
moduleTraitCor 
write.csv(moduleTraitCor,paste0(title,"-moduleTraitCor -MT.csv"))

save(design, textMatrix,moduleTraitCor, file = paste0(title,"-MT-heatmap.Rdata"))

################################################################################
# Set the title for the analysis
title = "GSE179285"

# Load the WGCNA data and module information
load("~/AD/IBD/result/step1_input.Rdata")
load("~/AD/IBD/result/step8_TOM.Rdata")
load("~/AD/IBD/result/step3_stepByStep_genes_modules.Rdata")

# Choose the first module from the moduleTraitCor for processing
choose_module = sub("ME","",rownames(moduleTraitCor)[1])

# Loop through the chosen module
for (i in choose_module){
  module = i
  cat("PROCESS:",title ,"\n")  # Print the processing title
  
  # Calculate the KME (module eigengene-based connectivity) values for the module
  datKME=signedKME(datExpr0, MEs, outputColumnName="kME_MM.")
  
  # Extract the module names from the MEs object
  modNames = substring(names(MEs), 3)
  
  # Find the column index for the module
  column = match(module, modNames)
  
  # Identify the genes that belong to the current module
  moduleGenes = moduleColors==module
  
  # Create a data frame of the genes in the module
  lightyellow_module<-as.data.frame(dimnames(data.frame(datExpr0))[[2]][moduleGenes])
  names(lightyellow_module)="genename"
  
  # Extract the KME values for the genes in the module
  lightyellow_KME<-as.data.frame(datKME[moduleGenes,column]) 
  names(lightyellow_KME)="KME"
  
  # Set the row names of the KME data frame to the gene names
  rownames(lightyellow_KME)=lightyellow_module$genename
  
  # Filter genes with KME values greater than 0.8
  FilterGenes = abs(lightyellow_KME$KME) > 0.8
  table(FilterGenes)  # Display the count of genes meeting the KME threshold
  
  # Subset the KME data frame to only include genes with KME values greater than 0.8
  lightyellow_hub<-subset(lightyellow_KME, abs(lightyellow_KME$KME)>0.8)
  
  # Write the hub genes to a CSV file
  write.csv(lightyellow_hub,paste0(title,paste("-hubgene_KME_",module,".csv")))
  
  # Print a message indicating the export of the hub genes
  cat("EXPORT:",module," hubgene","\n")
}
