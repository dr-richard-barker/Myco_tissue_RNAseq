#!/usr/bin/env Rscript
# Phase 4 -- normalisation and differential expression, per NASA GeneLab GL-DPPD-7101-G.
#
# Two parallel tracks, as the G revision of the pipeline specifies:
#   all-genes  : every feature
#   rRNA-removed: rRNA-biotype features dropped, then INDEPENDENTLY re-normalised
# The second is not cosmetic here. rRNA is 20-72% of reads and its share varies with tissue
# (Exudophore ~62% vs Nodule ~46%), so leaving it in imposes a composition bias that is
# confounded with the factor of interest.
#
# No TPM/FPKM anywhere: this is 3'-end tag counting, one count per molecule regardless of
# transcript length, so effective-length normalisation is meaningless.

suppressPackageStartupMessages({
  library(DESeq2); library(edgeR)
})

root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
args <- commandArgs(TRUE)
counts_file <- if (length(args) >= 1) args[1] else file.path(root, "counts", "counts_dedup.txt")
min_mrna    <- if (length(args) >= 2) as.numeric(args[2]) else 1e5
# drop_low=TRUE runs the sensitivity analysis: exclude samples below min_mrna entirely
# rather than only flagging them. The PCA shows samples under ~20k non-rRNA counts scatter
# while everything above clusters by tissue, so this separates signal from noise-dominated
# libraries instead of letting the latter dictate dispersion estimates.
drop_low    <- length(args) >= 3 && tolower(args[3]) %in% c("true", "drop", "1")
label       <- if (length(args) >= 4) args[4] else "PC9.15"
outdir      <- file.path(root, "results",
                         paste0("dge_", label, if (drop_low) "_filtered" else ""))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- load counts -------------------------------------------------------------
fc <- read.delim(counts_file, comment.char = "#", check.names = FALSE)
mat <- as.matrix(fc[, 7:ncol(fc)])
rownames(mat) <- fc$Geneid
colnames(mat) <- sub("\\..*$", "", basename(colnames(mat)))

# ---- biotypes from the final GTF --------------------------------------------
gtf <- file.path(root, "refs", label, paste0(label, "_final.gtf"))
gl <- readLines(gtf)
gl <- gl[!startsWith(gl, "#")]
f <- vapply(strsplit(gl, "\t"), `[`, "", 3)
gl <- gl[f == "gene"]
gid <- sub('.*gene_id "([^"]+)".*', "\\1", gl)
bio <- ifelse(grepl('gene_biotype "', gl), sub('.*gene_biotype "([^"]+)".*', "\\1", gl), "unknown")
biotype <- setNames(bio, gid)

# ---- metadata ---------------------------------------------------------------
rs <- read.csv(file.path(root, "metadata", "runsheet.csv"), check.names = FALSE)
rs$sample_name <- as.character(rs$sample_name)
stopifnot(all(colnames(mat) %in% rs$sample_name))
rs <- rs[match(colnames(mat), rs$sample_name), ]
coldata <- data.frame(
  row.names = colnames(mat),
  Tissue = factor(make.names(rs[["Factor Value[Tissue]"]])),
  well   = rs$well
)

is_rrna <- biotype[rownames(mat)] %in% c("rRNA")
mrna_totals <- colSums(mat[!is_rrna, , drop = FALSE])
cat("assigned non-rRNA counts per sample:\n"); print(mrna_totals)

low <- names(mrna_totals)[mrna_totals < min_mrna]
if (length(low)) {
  cat("\nFLAGGED as low-yield (<", format(min_mrna, big.mark = ","), " non-rRNA counts):\n", sep = "")
  print(data.frame(sample = low, counts = mrna_totals[low],
                   tissue = as.character(coldata[low, "Tissue"]), row.names = NULL))
}
write.csv(data.frame(sample = names(mrna_totals), mrna_counts = mrna_totals,
                     tissue = as.character(coldata$Tissue),
                     flagged_low = names(mrna_totals) %in% low, row.names = NULL),
          file.path(outdir, "sample_yield.csv"), row.names = FALSE)

if (drop_low && length(low)) {
  keep_s <- setdiff(colnames(mat), low)
  mat <- mat[, keep_s, drop = FALSE]
  coldata <- droplevels(coldata[keep_s, , drop = FALSE])
  cat("\nSENSITIVITY RUN: dropped", length(low), "samples; retained", ncol(mat), "\n")
  print(table(coldata$Tissue))
  if (any(table(coldata$Tissue) < 2)) {
    cat("WARNING: a tissue has fewer than 2 replicates; its contrasts are not interpretable\n")
  }
}

# ---- one track ---------------------------------------------------------------
run_track <- function(m, track) {
  cat("\n=====", track, "=====\n")
  keep <- filterByExpr(DGEList(m), group = coldata$Tissue)
  cat("genes retained by filterByExpr:", sum(keep), "of", nrow(m), "\n")
  m <- m[keep, , drop = FALSE]
  if (nrow(m) < 50) { cat("too few genes retained; skipping track\n"); return(invisible(NULL)) }

  dds <- DESeqDataSetFromMatrix(m, coldata, ~ Tissue)
  dds <- DESeq(dds, quiet = TRUE)

  write.csv(as.data.frame(counts(dds, normalized = TRUE)),
            file.path(outdir, paste0("normalized_counts_", track, ".csv")))
  vsd <- tryCatch(vst(dds, blind = TRUE),
                  error = function(e) varianceStabilizingTransformation(dds, blind = TRUE))
  write.csv(as.data.frame(assay(vsd)), file.path(outdir, paste0("VST_counts_", track, ".csv")))

  # PCA -- the plan's sanity gate is that replicates cluster by tissue
  pca <- prcomp(t(assay(vsd)))
  pv <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  pcs <- data.frame(sample = colnames(m), PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                    Tissue = coldata$Tissue, mrna = colSums(m))
  write.csv(pcs, file.path(outdir, paste0("PCA_", track, ".csv")), row.names = FALSE)
  cat("PC1", pv[1], "% PC2", pv[2], "%\n")
  cat("cor(PC1, log10 non-rRNA depth) =",
      round(cor(pcs$PC1, log10(pcs$mrna + 1)), 3), " <- depth artefact check\n")

  # LRT: any tissue effect at all
  lrt <- DESeq(dds, test = "LRT", reduced = ~ 1, quiet = TRUE)
  res_lrt <- results(lrt)
  cat("LRT any-tissue effect, padj<0.05:", sum(res_lrt$padj < 0.05, na.rm = TRUE), "genes\n")
  write.csv(as.data.frame(res_lrt[order(res_lrt$padj), ]),
            file.path(outdir, paste0("LRT_", track, ".csv")))

  lv <- levels(coldata$Tissue)
  for (i in seq_along(lv)) for (j in seq_len(i - 1)) {
    r <- results(dds, contrast = c("Tissue", lv[i], lv[j]))
    n <- sum(r$padj < 0.05 & abs(r$log2FoldChange) > 1, na.rm = TRUE)
    cat(sprintf("  %-22s vs %-22s : %5d DE (padj<0.05, |LFC|>1)\n", lv[i], lv[j], n))
    write.csv(as.data.frame(r[order(r$padj), ]),
              file.path(outdir, sprintf("DE_%s_%s_vs_%s.csv", track, lv[i], lv[j])))
  }
  invisible(dds)
}

run_track(mat, "all_genes")
run_track(mat[!is_rrna, , drop = FALSE], "rRNArm")

cat("\nresults written to", outdir, "\n")
