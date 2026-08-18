#!/usr/bin/env Rscript
# Phase 5 -- weighted gene co-expression network analysis (WGCNA).
#
# POWER CAVEAT, stated up front. WGCNA's own guidance is that fewer than ~15 samples gives
# unstable modules; this network has 12, across 4 groups with n=2-4. Modules here are
# EXPLORATORY. They are reported with module sizes, soft-threshold diagnostics and
# module-trait correlations so the reader can judge them, and hub genes are cross-referenced
# to CAZyme/secretome annotation rather than interpreted on their own.
#
# Input is the all_genes VST (4,752 genes) with rRNA features removed, rather than the
# rRNArm track (1,666): WGCNA behaves better with more genes, and the rRNA features are
# dropped explicitly here anyway.

suppressPackageStartupMessages({ library(WGCNA) })
options(stringsAsFactors = FALSE)
enableWGCNAThreads(8)

root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
label <- "BOM_ss5"
outdir <- file.path(root, "results", "wgcna"); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

vst <- read.csv(file.path(root, "results", paste0("dge_", label, "_filtered"), "VST_counts_all_genes.csv"), row.names = 1, check.names = FALSE)

# drop rRNA features by biotype
gtf <- readLines(file.path(root, "refs", label, paste0(label, "_final.gtf")))
gtf <- gtf[!startsWith(gtf, "#")]
f3 <- vapply(strsplit(gtf, "\t"), `[`, "", 3); gtf <- gtf[f3 == "gene"]
gid <- sub('.*gene_id "([^"]+)".*', "\\1", gtf)
bio <- ifelse(grepl('gene_biotype "', gtf), sub('.*gene_biotype "([^"]+)".*', "\\1", gtf), "unknown")
rrna <- gid[bio == "rRNA"]
vst <- vst[!rownames(vst) %in% rrna, , drop = FALSE]
cat("genes after rRNA removal:", nrow(vst), " samples:", ncol(vst), "\n")

datExpr <- as.data.frame(t(vst))
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  cat("dropping", sum(!gsg$goodGenes), "genes and", sum(!gsg$goodSamples), "samples failing QC\n")
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}
cat("final matrix:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n")

# ---- soft threshold ----
powers <- c(1:10, seq(12, 30, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)
write.csv(sft$fitIndices, file.path(outdir, "soft_threshold.csv"), row.names = FALSE)
best <- sft$powerEstimate
if (is.na(best)) {
  # WGCNA's default when no power reaches R^2 0.9 -- for signed networks with <20 samples
  # the authors recommend 18; take that rather than silently using a bad fit
  best <- 18
  cat("NOTE: no power reached scale-free R^2 0.9 (expected at n=12); using recommended", best, "\n")
}
cat("soft-threshold power:", best, "\n")
r2 <- sft$fitIndices[sft$fitIndices$Power == best, "SFT.R.sq"]
cat("scale-free topology R^2 at chosen power:", round(r2, 3), "\n")

# ---- modules ----
net <- blockwiseModules(datExpr, power = best, networkType = "signed",
                        TOMType = "signed", minModuleSize = 30,
                        mergeCutHeight = 0.25, numericLabels = TRUE,
                        maxBlockSize = 6000, verbose = 0)
moduleColors <- labels2colors(net$colors)
tab <- table(moduleColors)
cat("\nmodules found:", length(tab) - ("grey" %in% names(tab)), "(+ grey = unassigned)\n")
print(sort(tab, decreasing = TRUE))

write.csv(data.frame(gene = colnames(datExpr), module = moduleColors,
                     row.names = NULL),
          file.path(outdir, "gene_modules.csv"), row.names = FALSE)

# ---- module-trait correlation ----
rs <- read.csv(file.path(root, "metadata", "runsheet.csv"), check.names = FALSE)
tis <- rs[["Factor Value[Tissue]"]][match(rownames(datExpr), rs$sample_name)]
traits <- model.matrix(~ 0 + factor(tis)); colnames(traits) <- levels(factor(tis))

MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
modTraitCor <- cor(MEs, traits, use = "p")
modTraitP <- corPvalueStudent(modTraitCor, nrow(datExpr))
write.csv(modTraitCor, file.path(outdir, "module_trait_correlation.csv"))
write.csv(modTraitP, file.path(outdir, "module_trait_pvalue.csv"))

cat("\nmodule-trait correlations (r, p) -- only |r|>0.7 & p<0.05 shown:\n")
hits <- 0
for (m in rownames(modTraitCor)) for (t in colnames(modTraitCor)) {
  if (abs(modTraitCor[m, t]) > 0.7 && modTraitP[m, t] < 0.05) {
    cat(sprintf("  %-22s %-20s r=%+.2f  p=%.4f  (n genes=%d)\n", m, t,
                modTraitCor[m, t], modTraitP[m, t],
                sum(moduleColors == sub("^ME", "", m))))
    hits <- hits + 1
  }
}
if (hits == 0) cat("  none -- no module is confidently tissue-associated at this sample size\n")

# ---- hub genes: highest module membership within each module ----
MM <- as.data.frame(cor(datExpr, MEs, use = "p"))
hub <- do.call(rbind, lapply(setdiff(unique(moduleColors), "grey"), function(m) {
  g <- colnames(datExpr)[moduleColors == m]
  k <- MM[g, paste0("ME", m)]
  o <- order(-k)[1:min(20, length(g))]
  data.frame(module = m, gene = g[o], module_membership = round(k[o], 3))
}))
write.csv(hub, file.path(outdir, "hub_genes.csv"), row.names = FALSE)
cat("\nwrote", outdir, "\n")
