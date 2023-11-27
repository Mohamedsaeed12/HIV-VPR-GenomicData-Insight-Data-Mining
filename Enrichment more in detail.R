egoMF <- enrichGO(gene         = significant_gene_symbols,
                OrgDb        = org.Hs.eg.db,
                keyType      = "SYMBOL",
                ont          = "MF", # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                pAdjustMethod = "BH", # Molecular function
                qvalueCutoff = 0.05,
                readable     = TRUE)

egoBP <- enrichGO(gene         = significant_gene_symbols,
                OrgDb        = org.Hs.eg.db,
                keyType      = "SYMBOL",
                ont          = "BP", # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                pAdjustMethod = "BH", # Molecular function
                qvalueCutoff = 0.05,
                readable     = TRUE)

egoCC <- enrichGO(gene         = significant_gene_symbols,
                OrgDb        = org.Hs.eg.db,
                keyType      = "SYMBOL",
                ont          = "CC", # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                pAdjustMethod = "BH", # Molecular function
                qvalueCutoff = 0.05,
                readable     = TRUE)



dotplot(ego)
dotplot(egoMF)
dotplot(egoBP)

#extracting results from egoBP
egoBP_results <- egoBP@result
head(egoBP_results)

library(ggplot2)

ggplot(egoBP_results, aes(x = reorder(Description, -pvalue), y = -log10(pvalue))) +
  geom_point() +
  coord_flip() +
  labs(x = "GO Term", y = "-log10(p-value)", title = "GO Enrichment Analysis") +
  theme_minimal()

#Since the plot is clustered we need to do some filtering to find the top biological processes
#P-value filteration

# Adjust the threshold as needed
pvalue_threshold = 0.05

filtered_egoBP_results <- egoBP_results[egoBP_results$pvalue < pvalue_threshold, ]

# Then plot using the filtered data
ggplot(filtered_egoBP_results, aes(x = reorder(Description, -pvalue), y = -log10(pvalue))) +
  geom_point() +
  coord_flip() +
  labs(x = "GO Term", y = "-log10(p-value)", title = "GO Enrichment Analysis") +
  theme_minimal()


#N Term filtration

top_n = 20  # Adjust the number as needed

top_terms_egoBP_results <- egoBP_results[order(egoBP_results$pvalue), ][1:top_n, ]

# Plot using the top terms
ggplot(top_terms_egoBP_results, aes(x = reorder(Description, -pvalue), y = -log10(pvalue))) +
  geom_point() +
  coord_flip() +
  labs(x = "Gene ontology", y = "-log10(p-value)", title = "Top 20 GO Enrichment Terms") +
  theme_minimal()


#extracting results from egoBP
egoCC_results <- egoCC@result


library(ggplot2)

ggplot(egoCC_results, aes(x = reorder(Description, -pvalue), y = -log10(pvalue))) +
  geom_point() +
  coord_flip() +
  labs(x = "GO Term", y = "-log10(p-value)", title = "GO Enrichment Analysis") +
  theme_minimal()

#Since the plot is clustered we need to do some filtering to find the top biological processes
#P-value filteration

# Adjust the threshold as needed
pvalue_threshold = 0.05

filtered_egoCC_results <- egoCC_results[egoCC_results$pvalue < pvalue_threshold, ]

# Then plot using the filtered data
ggplot(filtered_egoCC_results, aes(x = reorder(Description, -pvalue), y = -log10(pvalue))) +
  geom_point() +
  coord_flip() +
  labs(x = "GO Term", y = "-log10(p-value)", title = "GO Enrichment Analysis") +
  theme_minimal()


#N Term filtration

top_n = 20  # Adjust the number as needed

top_terms_egoCC_results <- egoCC_results[order(egoCC_results$pvalue), ][1:top_n, ]

# Plot using the top terms
ggplot(top_terms_egoCC_results, aes(x = reorder(Description, -pvalue), y = -log10(pvalue))) +
  geom_point() +
  coord_flip() +
  labs(x = "Gene ontology", y = "-log10(p-value)", title = "Top 20 GO Enrichment Terms") +
  theme_minimal()


