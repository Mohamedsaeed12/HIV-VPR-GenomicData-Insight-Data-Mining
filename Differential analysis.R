# Differential analysis

library(limma)

# Correcting the conditions based on your data structure
conditions <- rep(c("Vpr-Mac", "Zs-Mac"), ncol(expr_data) / 2)
design <- model.matrix(~0 + conditions)
colnames(design) <- levels(conditions)

fit <- lmFit(expr_data, design)

print(design)

colnames(design) <- c("VprMac", "ZsMac")

contrast.matrix <- makeContrasts(Difference = ZsMac - VprMac, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, coef="Difference", adjust="fdr", number=Inf)

head(results)


#logFC: This is the log2 fold change. A positive value indicates that the gene is upregulated in "ZsMac" compared to "VprMac", while a negative value indicates downregulation.

#AveExpr: This is the average expression of the gene across all samples.

#t: This is the t-statistic value, which indicates the difference in expression of a gene between the two conditions relative to the variability in expression across samples.

#P.Value: This is the p-value for the hypothesis test that checks if the gene is differentially expressed between the two conditions. A small p-value (typically < 0.05) suggests that the gene is differentially expressed.

#adj.P.Val: This is the adjusted p-value, which corrects for multiple testing using the Benjamini & Hochberg method. It controls the false discovery rate (FDR). Genes with an adjusted p-value below a certain threshold (e.g., 0.05) are considered to be differentially expressed.

#B: This is the B-statistic or log-odds that the gene is differentially expressed. A higher B value indicates stronger evidence for differential expression.

#From the head(results) output, you can see the top differentially expressed genes ranked by their t-statistic. The genes at the top of this list show the most significant differences in expression between the two conditions.