#!/usr/bin/env Rscript
# Produce DESeq2 size-factor normalised counts for EVERY gene, with no expression filter.
#
# The differential-expression matrices are filtered with edgeR::filterByExpr, which is correct
# for hypothesis testing but wrong as input to metabolic contextualisation: a gene absent from
# the matrix is indistinguishable from a gene expressed at zero, and RIPTiDe then prunes every
# reaction that depends on it. Using the filtered matrix pruned the models from 5,247 reactions
# to ~90. This writes the full matrix (rRNA features removed, 12 retained libraries).

suppressPackageStartupMessages(library(DESeq2))
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
label <- "BOM_ss5"

fc <- read.delim(file.path(root, "counts", paste0("counts_", label, "_dedup.txt")),
                 comment.char = "#", check.names = FALSE)
mat <- as.matrix(fc[, 7:ncol(fc)]); rownames(mat) <- fc$Geneid
colnames(mat) <- sub("\\..*$", "", basename(colnames(mat)))

gtf <- readLines(file.path(root, "refs", label, paste0(label, "_final.gtf")))
gtf <- gtf[!startsWith(gtf, "#")]
f3 <- vapply(strsplit(gtf, "\t"), `[`, "", 3); gtf <- gtf[f3 == "gene"]
gid <- sub('.*gene_id "([^"]+)".*', "\\1", gtf)
bio <- ifelse(grepl('gene_biotype "', gtf), sub('.*gene_biotype "([^"]+)".*', "\\1", gtf), "unknown")
rrna <- gid[bio == "rRNA"]
mat <- mat[!rownames(mat) %in% rrna, , drop = FALSE]

yield <- read.csv(file.path(root, "results", paste0("dge_", label, "_filtered"), "sample_yield.csv"))
keep <- yield$sample[!yield$flagged_low]
mat <- mat[, colnames(mat) %in% keep, drop = FALSE]
cat("matrix:", nrow(mat), "genes x", ncol(mat), "libraries\n")

rs <- read.csv(file.path(root, "metadata", "runsheet.csv"), check.names = FALSE)
cd <- data.frame(row.names = colnames(mat),
                 Tissue = factor(make.names(rs[["Factor Value[Tissue]"]][match(colnames(mat), rs$sample_name)])))
dds <- DESeqDataSetFromMatrix(mat, cd, ~ Tissue)
dds <- estimateSizeFactors(dds)
norm <- counts(dds, normalized = TRUE)
out <- file.path(root, "results", "normalised_all_genes.csv")
write.csv(as.data.frame(norm), out)
cat("wrote", out, "\n")
cat("genes with any expression:", sum(rowSums(norm) > 0), "\n")
