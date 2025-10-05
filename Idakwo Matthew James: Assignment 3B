if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
pkgs <- c("affy", "affyPLM", "arrayQualityMetrics", "limma", "genefilter", 
          "annotate", "hgu133plus2.db", "GEOquery", "ggplot2", "pheatmap", "R.utils")
for(p in pkgs) if(!requireNamespace(p, quietly=TRUE)) BiocManager::install(p)

library(GEOquery)
library(affy)
library(affyPLM)
library(arrayQualityMetrics)
library(limma)
library(genefilter)
library(annotate)
library(hgu133plus2.db)
library(ggplot2)
library(pheatmap)
library(R.utils)

setwd(tempdir())
getGEOSuppFiles("GSE42394")

untar("GSE42394/GSE42394_RAW.tar", exdir = "GSE42394/CEL")
celfiles <- list.files("GSE42394/CEL", pattern = "\\.gz$", full.names = TRUE)
sapply(celfiles, R.utils::gunzip, overwrite = TRUE)

raw <- ReadAffy(celfile.path = "GSE42394/CEL")

gse <- getGEO("GSE42394", GSEMatrix = TRUE)
pheno <- pData(gse[[1]])
rownames(pheno) <- sampleNames(raw)
pheno$group <- ifelse(grepl("Glioblastoma", pheno$title, ignore.case = TRUE), 
                      "cancer", "normal")

exprs_raw <- exprs(raw)
pdf("QC_boxplot_raw.pdf")
boxplot(raw, main = "Raw Intensities - Before Normalization", las = 2)
dev.off()

pca_raw <- prcomp(t(exprs_raw), scale. = TRUE)
pdf("QC_PCA_raw.pdf")
plot(pca_raw$x[,1], pca_raw$x[,2], col = as.factor(pheno$group),
     main = "PCA Before Normalization", pch = 19)
legend("topright", legend = unique(pheno$group), col = 1:2, pch = 19)
dev.off()

med_raw <- apply(exprs_raw, 2, median)
z_raw <- scale(med_raw)
outliers_before <- names(which(abs(z_raw) > 2))
cat("Outliers before normalization:", outliers_before, "\n")
cat("Count:", length(outliers_before), "\n")

eset <- rma(raw)
exprs_norm <- exprs(eset)

pdf("QC_boxplot_normalized.pdf")
boxplot(exprs_norm, main = "Normalized Intensities - After RMA", las = 2)
dev.off()

pca_norm <- prcomp(t(exprs_norm), scale. = TRUE)
pdf("QC_PCA_normalized.pdf")
plot(pca_norm$x[,1], pca_norm$x[,2], col = as.factor(pheno$group),
     main = "PCA After Normalization", pch = 19)
legend("topright", legend = unique(pheno$group), col = 1:2, pch = 19)
dev.off()

med_norm <- apply(exprs_norm, 2, median)
z_norm <- scale(med_norm)
outliers_after <- names(which(abs(z_norm) > 2))
cat("Outliers after normalization:", outliers_after, "\n")
cat("Count:", length(outliers_after), "\n")

n_before <- nrow(exprs_norm)
keep <- rowMeans(exprs_norm) > log2(100)
exprs_filtered <- exprs_norm[keep, ]
n_after <- nrow(exprs_filtered)
cat("Number of probes before filtering:", n_before, "\n")
cat("Number of probes after filtering:", n_after, "\n")

cat("\nSUMMARY:\n")
cat("Arrays flagged BEFORE normalization:", length(outliers_before), "\n")
cat("Arrays flagged AFTER normalization:", length(outliers_after), "\n")
cat("Probes before filtering:", n_before, "\n")
cat("Probes after filtering:", n_after, "\n")

write.table(exprs_filtered, "GSE42394_filtered_expression.txt", sep="\t", quote=FALSE)
write.csv(pheno, "GSE42394_pheno.csv", row.names=TRUE)
