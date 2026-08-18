#!/usr/bin/env Rscript
# Phase 5 -- is the exudophore "SSU-rRNA maturation" module real biology or an rRNA artefact?
#
# The lightcyan module correlates with Exudophore at r=+0.91 and enriches for SSU-rRNA
# maturation. These libraries are 20-70% rRNA, and rRNA fraction varies systematically by
# tissue, so a module of ribosome-biogenesis genes could reflect either genuine high growth or
# co-variation with whatever technical factor drives rRNA content.
#
# Test: for every module eigengene, compare
#   (a) correlation with the per-sample rRNA fraction  -- the technical axis
#   (b) correlation with tissue membership             -- the biological axis
#   (c) PARTIAL correlation with tissue, controlling for rRNA fraction
# A module whose tissue association survives (c) is not explained by rRNA load.

suppressPackageStartupMessages({ library(WGCNA) })
options(stringsAsFactors = FALSE)

root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
label <- "BOM_ss5"

vst <- read.csv(file.path(root, "results", paste0("dge_", label, "_filtered"), "VST_counts_all_genes.csv"),
                row.names = 1, check.names = FALSE)
mods <- read.csv(file.path(root, "results/wgcna/gene_modules.csv"))
budget <- read.csv(file.path(root, "results/read_budget_BOM_ss5.csv"))
rs <- read.csv(file.path(root, "metadata", "runsheet.csv"), check.names = FALSE)

vst <- vst[rownames(vst) %in% mods$gene, , drop = FALSE]
datExpr <- as.data.frame(t(vst))
mc <- mods$module[match(colnames(datExpr), mods$gene)]

MEs <- orderMEs(moduleEigengenes(datExpr, mc)$eigengenes)
samples <- rownames(datExpr)

# technical axis: rRNA as a fraction of assigned counts
b <- budget[match(samples, budget$sample), ]
rrna_frac <- b$rRNA / (b$rRNA + b$mRNA)
depth <- log10(b$mRNA + 1)
tis <- rs[["Factor Value[Tissue]"]][match(samples, rs$sample_name)]

cat("samples:", length(samples), "\n")
cat("rRNA fraction range:", sprintf("%.3f-%.3f", min(rrna_frac), max(rrna_frac)), "\n")
cat("cor(rRNA fraction, log10 mRNA depth) =", round(cor(rrna_frac, depth), 3), "\n\n")

partial <- function(x, y, z) {
  rxy <- cor(x, y); rxz <- cor(x, z); ryz <- cor(y, z)
  (rxy - rxz * ryz) / sqrt((1 - rxz^2) * (1 - ryz^2))
}
pval <- function(r, n, k = 0) {
  df <- n - 2 - k
  t <- r * sqrt(df / (1 - r^2))
  2 * pt(-abs(t), df)
}

out <- list()
for (t in unique(tis)) {
  trait <- as.numeric(tis == t)
  for (m in colnames(MEs)) {
    e <- MEs[[m]]
    r_tis  <- cor(e, trait)
    r_rrna <- cor(e, rrna_frac)
    r_part <- partial(e, trait, rrna_frac)
    out[[length(out) + 1]] <- data.frame(
      module = sub("^ME", "", m), tissue = t,
      r_tissue = r_tis, p_tissue = pval(r_tis, length(samples)),
      r_rRNAfrac = r_rrna,
      r_tissue_partial = r_part, p_partial = pval(r_part, length(samples), 1))
  }
}
res <- do.call(rbind, out)
res$fdr_tissue  <- p.adjust(res$p_tissue,  "BH")
res$fdr_partial <- p.adjust(res$p_partial, "BH")
write.csv(res, file.path(root, "results/wgcna/rrna_confound_test.csv"), row.names = FALSE)

sig <- res[res$fdr_tissue < 0.05 & abs(res$r_tissue) > 0.7, ]
sig <- sig[order(sig$fdr_tissue), ]
cat("modules with a tissue association surviving FDR (from the original test):\n\n")
cat(sprintf("%-18s %-20s %8s %10s %10s %10s %s\n",
            "module", "tissue", "r_tissue", "r_rRNAfrac", "r_partial", "FDR_part", "verdict"))
for (i in seq_len(nrow(sig))) {
  s <- sig[i, ]
  verdict <- if (s$fdr_partial < 0.05 && abs(s$r_tissue_partial) > 0.7) "survives" else
             if (abs(s$r_rRNAfrac) > 0.7) "CONFOUNDED by rRNA" else "weakened"
  cat(sprintf("%-18s %-20s %+8.2f %+10.2f %+10.2f %10.4f %s\n",
              s$module, s$tissue, s$r_tissue, s$r_rRNAfrac,
              s$r_tissue_partial, s$fdr_partial, verdict))
}
cat("\nwrote results/wgcna/rrna_confound_test.csv\n")
