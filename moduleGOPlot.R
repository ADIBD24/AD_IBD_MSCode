####
# Read gene names from a file
genes <- nodeFile$nodeName

# Load the node file which contains gene names
nodeFile <- read.csv("step8_CytoscapeInput-nodes-lightpink4.txt", sep = "\t")

# Perform GO enrichment analysis
enrich.go <- enrichGO(
  gene = genes,  # Gene list for enrichment
  OrgDb = org.Hs.eg.db,  # Organism database (here, sheep is used as an example)
  keyType = 'SYMBOL',  # Gene name type (here, it's assumed to be Entrez IDs)
  ont = 'ALL',  # GO Ontology categories (BP, MF, CC) or 'ALL' for all
  pAdjustMethod = "BH",  # P-value adjustment method
  pvalueCutoff = 1,  # P-value cutoff (set to 1 to output all)
  qvalueCutoff = 1,  # Q-value cutoff (set to 1 to output all)
  readable = FALSE
)

# Save the GO enrichment results to a CSV file
write.csv(enrich.go, paste(paste0(title,"-GO-"), paste(module, collapse="-"), ".csv", sep=""), row.names = FALSE, quote = FALSE)
# Convert the enrichment result to a data frame
enrich.go <- as.data.frame(enrich.go)

# Format the GO terms and create a factor for ordered plotting
enrich.go$Term <- paste(enrich.go$ID, enrich.go$Description, sep = ': ')
enrich.go$Term <- factor(enrich.go$Term, levels = enrich.go$Term, ordered = T)

# Create a vertical bar plot for GO enrichment
p1 <- ggplot(enrich.go,
             aes(x=Term, y=Count, fill=ONTOLOGY)) +
  geom_bar(stat="identity", width=0.8) +
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666")) +
  coord_flip() +
  xlab("GO term") +
  ylab("Gene Number") +
  labs(title = "GO Terms Enrich") +
  theme_bw()
# Add facet grid based on ONTOLOGY categories
p1 + facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')

# Save the plot as a PDF file
ggsave(paste(paste0(title,"-GO-"), paste(module, collapse="-"), ".pdf", sep=""))

# Create a horizontal bar plot for GO enrichment based on ONTOLOGY types
p2 <- ggplot(eGo10,
             aes(x=Term, y=Count, fill=ONTOLOGY)) +
  geom_bar(stat="identity", width=0.8) +
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666")) +
  xlab("GO term") +
  ylab("Gene Number") +
  labs(title = "GO Terms Enrich") +
  theme_bw() +
  theme(axis.text.x=element_text(family="sans", face = "bold", color="gray50", angle = 70, vjust = 1, hjust = 1))
# Add facet grid based on ONTOLOGY categories
p2 + facet_grid(.~ONTOLOGY, scale = 'free_x', space = 'free_x')

# Save the horizontal bar plot as a PDF file with specified dimensions
ggsave(paste(paste0(title,"-GO-"), paste(module, collapse="-"), ".pdf", sep=""), width = 20, height = 10)

# Load tidyverse and ggplot2 libraries for further data manipulation and plotting
library(tidyverse)
library(ggplot2)

# Convert the enrichment results to a data frame and manipulate it
eGo <- as.data.frame(enrich.go)
# Split the GeneRatio and BgRatio columns into two separate columns each
eGo <- separate(data=eGo, col=GeneRatio, into = c("GR1", "GR2"), sep = "/")
eGo <- separate(data=eGo, col=BgRatio, into = c("BR1", "BR2"), sep = "/")
# Calculate the Enrichment Factor
eGo <- mutate(eGo, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2)))

# Filter top 10 GO terms for each ontology (Biological Process, Cellular Component, Molecular Function)
eGoBP <- eGo %>% filter(ONTOLOGY=="BP") %>% filter(row_number() >= 1, row_number() <= 10)
eGoCC <- eGo %>% filter(ONTOLOGY=="CC") %>% filter(row_number() >= 1, row_number() <= 10)
eGoMF <- eGo %>% filter(ONTOLOGY=="MF") %>% filter(row_number() >= 1, row_number() <= 10)
eGo10 <- rbind(eGoBP, eGoMF, eGoCC)

# Create a scatter plot for the top 10 GO terms of each ontology
p <- ggplot(eGo10, aes(enrichment_factor, fct_reorder(factor(Description), enrichment_factor))) +
  geom_point(aes(size=Count, color=-1*log10(pvalue), shape=ONTOLOGY)) +
  scale_color_gradient(low="#5ab4ac", high ="#fc8d59") +
  labs(color=expression(-log[10](p-value)), size="Count", shape="Ontology",
       x="Enrichment Factor", y="GO term", title="GO enrichment") +
  theme_bw() + 
  facet_wrap( ~ ONTOLOGY, ncol= 1, scale='free')
# Save the scatter plot as a PDF file with specified dimensions
ggsave(paste(paste0(title,"-GO-"), paste(module, collapse="-"), ".pdf", sep=""), width = 8, height = 10)

# Load GOplot library for GO term visualization
library(GOplot)
# Rename columns of the GO enrichment results for compatibility with GOplot functions
colnames(enrich.go) <- c("Category", "ID", "Term", "GeneRatio", "BgRatio", "pvalue", "adj_pval", "qvalue", "Genes", "Count")
# Replace the "/" separator in the Genes column with a comma
enrich.go$Genes <- gsub("/", ",", enrich.go$Genes)

# Prepare data for GOplot visualizations
GOterms <- data.frame(
  category = eGo10$ONTOLOGY,
  ID = eGo10$ID,
  term = eGo10$Description, 
  genes = gsub("/", ",", eGo10$geneID), 
  adj_pval = eGo10$p.adjust
)

# Load DEG data and prepare it for GOplot
load("~/AD/The RNAseq Harmonization Study/result/DEG/MSBB_CQN-deg.Rdata")
colnames(DEG)[1] <- "ID"
circ <- circle_dat(GOterms, as.data.frame(DEG))

# Create a bubble plot for GO terms and save it as a PDF file
GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)
ggsave(paste(paste0(title,"-GOBubble-"), paste(module, collapse="-"), ".pdf", sep=""))

# Reduce redundant terms and create a bubble plot with a label threshold
reduced_circ <- reduce_overlap(circ, overlap = 0.75)
GOBubble(reduced_circ, labels = 2.8, display = 'multiple')
ggsave(paste(paste0(title,"-GOBubble-"), paste(module, collapse="-"), ".pdf", sep=""), width = 10, height = 8)

# Generate a circular plot for GO terms and save it as a PDF file
GOCircle(circ)
ggsave(paste(paste0(title,"-GOCircle-"), paste(module, collapse="-"), ".pdf", sep=""))

# Customize the circular plot with various aesthetic adjustments
# Save the customized circular plot as a PDF file with high resolution
ggsave(paste(paste0(title,"-GOCircle-"), paste(module, collapse="-"), ".pdf", sep=""), 
       dpi = 320, width = 14, height = 10)

# Create the plot
chord <- chord_dat(data = circ, genes = DEG, process = enrich.go$Term[1:8])
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)
ggsave(paste(paste0(title,"-GOChord-"), paste(module, collapse="-"), ".pdf", sep=""),width = 10,height = 8)

GOCluster(circ, enrich.go$Term[1:8], clust.by = 'term', lfc.col = c('darkgoldenrod1', 'black', 'cyan1'))
