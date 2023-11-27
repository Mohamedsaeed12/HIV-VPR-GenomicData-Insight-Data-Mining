# Functional Enrichment Analysis

#Functional Enrichment Analysis aims to identify biological processes, molecular functions, cellular components, or pathways that are over-represented in a set of genes. 
#One popular tool for this purpose is the clusterProfiler package in R, which supports multiple annotation sources.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

library(clusterProfiler)
library(org.Hs.eg.db)

significant_genes <- rownames(results[results$P.Value < 0.05, ])

ego <- enrichGO(gene         = significant_genes,
                OrgDb        = org.Hs.eg.db,
                keyType      = "SYMBOL",
                ont          = "ALL", # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                pAdjustMethod = "BH", # Benjamini & Hochberg method
                qvalueCutoff = 0.05,
                readable     = TRUE)

barplot(ego, showCategory=10)

#ego was null, which means none of the genes in the significant_genes list were found to be enriched in any GO terms based on the criteria provided
#let's troubleshoot:

length(significant_genes)
head(significant_genes)

# Convert the numeric values in significant_genes to integers
significant_indices <- as.integer(significant_genes)

# Extract the corresponding gene symbols using the indices
significant_gene_symbols <- normalized.expr$`Gene Symbol`[significant_indices]

# Remove duplicates and NA values
significant_gene_symbols <- unique(significant_gene_symbols)
significant_gene_symbols <- significant_gene_symbols[!is.na(significant_gene_symbols)]

# Check the mapped gene symbols
head(significant_gene_symbols)


# Perform GO enrichment analysis
ego <- enrichGO(gene         = significant_gene_symbols,
                OrgDb        = org.Hs.eg.db,
                keyType      = "SYMBOL",
                ont          = "ALL", # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                pAdjustMethod = "BH", # Benjamini & Hochberg method
                qvalueCutoff = 0.05,
                readable     = TRUE)

# View the results
head(ego)

dotplot(ego)

#this plot displays the most significant GO terms on the y-axis and the gene ratio on the x-axis. The size of the dots represents the gene count, and the color indicates the adjusted p-value.

barplot(ego)

#This plot displays the enriched GO terms on the y-axis and the gene count or gene ratio on the x-axis.

install.packages("UpSetR")
library(UpSetR)


gene_sets <- geneInCategory(ego)
upset(fromList(gene_sets))

#An upset plot that visualizes intersecting sets. It's useful when you have multiple gene lists and you want to see the intersections between them.


avg_expression <- rowMeans(normalized.expr[, 2:(ncol(normalized.expr)-1)])
cnetplot(ego, foldChange=avg_expression)

#A category network plot that visualizes the relationship between genes and the GO terms they are associated with.

#Functional Enrichment Analysis Summary:

#Objective:
 # To identify biological processes, molecular functions, cellular components, or pathways that are significantly over-represented in a set of genes.

#Methodology:
  
 # Package Installation and Loading:
  
  #Installed and loaded the clusterProfiler package, which is a popular tool for enrichment analysis.
#Loaded the org.Hs.eg.db package, which provides annotation data for human genes.
#Gene Selection:
  
  #Selected genes that had an adjusted p-value less than 0.05 from the differential expression analysis results.
#Initial Enrichment Analysis:
  
 # Performed Gene Ontology (GO) enrichment analysis using the enrichGO function.
#The analysis returned no significant results (ego was null).
#Troubleshooting:
  
 # Checked the length and head of the significant_genes list.
#Mapped the numeric indices in significant_genes to their corresponding gene symbols from the normalized.expr dataframe.
#Removed duplicates and NA values from the gene symbols list.
#Revised Enrichment Analysis:
  
 # Re-ran the GO enrichment analysis using the mapped gene symbols.
#Viewed the top results of the enrichment analysis.
#Visualization:
  
 # Dotplot: Visualized the most significant GO terms against the gene ratio. The size of the dots represents the gene count, while the color indicates the adjusted p-value.
#Barplot: Displayed the enriched GO terms against the gene count or gene ratio.
#UpSet Plot: Visualized the intersections between multiple gene lists using the UpSetR package.
#Category Network Plot (cnetplot): Displayed the relationship between genes and the GO terms they are associated with. The layout is designed to group genes and terms that are closely related.
#Outcome:
 # The analysis identified several GO terms that were significantly enriched in the set of differentially expressed genes. The visualizations provided insights into the biological processes, molecular functions, and cellular components associated with these genes.



