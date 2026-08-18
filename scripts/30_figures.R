#!/usr/bin/env Rscript
# Generate every manuscript figure as vector PDF at npj column widths.
#
# Sizing: npj/Nature Portfolio single column 88 mm, double column 180 mm. Output is PDF
# because the Springer Nature class is compiled here with XeTeX, which cannot embed EPS.
# Every panel reads a committed file under results/, so figures regenerate from a clean
# checkout without re-running any alignment.

suppressPackageStartupMessages({
  library(jsonlite); library(ggplot2); library(dplyr); library(tidyr); library(readr); library(patchwork); library(scales)
})
options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
FIG <- file.path(root, "latex", "figures"); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
MM <- function(x) x / 25.4          # mm -> inches
W1 <- MM(88); W2 <- MM(180)

# Okabe-Ito, colourblind safe
TIS <- c("Exuding mycelium" = "#0072B2", "Fuzzy mycelium" = "#009E73",
         "Exudophore"       = "#D55E00", "Nodule"         = "#CC79A7")
FATE <- c("rRNA" = "#B2182B", "mRNA" = "#2166AC", "multi-mapped" = "#F4A582",
          "no feature" = "#D1E5F0", "unmapped" = "#BDBDBD", "lost in QC" = "#EEEEEE")

theme_npj <- function(base = 7) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.2, colour = "grey92"),
          panel.border = element_rect(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.25),
          strip.background = element_rect(fill = "grey95", linewidth = 0.3),
          strip.text = element_text(size = base - 0.5, face = "bold"),
          legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0),
          plot.title = element_text(size = base + 0.5, face = "bold"),
          plot.tag = element_text(size = base + 2, face = "bold"))
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h,
         units = "in", device = cairo_pdf, limitsize = FALSE)
  cat(sprintf("  %-28s %.0f x %.0f mm\n", paste0(name, ".pdf"), w * 25.4, h * 25.4))
}

budget <- read_csv(file.path(root, "results/read_budget_BOM_ss5.csv"))
budget$tissue <- factor(budget$tissue, levels = names(TIS))
budget <- budget %>% arrange(tissue, well)
budget$well <- factor(budget$well, levels = budget$well)

# ------------------------------------------------------------------ Figure 1
# Read fate is reconciled exactly against raw reads: trimming + unaligned + UMI duplicates +
# featureCounts categories sum to the raw total (verified per library). An earlier version
# lumped trimming and duplicates into "lost in QC", which put 40-65% of reads in a category
# that is mostly successful deduplication rather than loss.
fig1 <- function() {
  smry <- file.path(root, "counts/counts_BOM_ss5_dedup.txt.summary")
  s <- read.delim(smry, check.names = FALSE); rownames(s) <- s$Status
  nm <- sub("\\..*$", "", basename(colnames(s)[-1]))
  get <- function(k) setNames(as.numeric(s[k, -1]), nm)
  dup <- read_csv(file.path(root, "results/figure_data/dedup_rates.csv"))

  trim <- sapply(budget$sample, function(x) {
    j <- jsonlite::fromJSON(file.path(root, "qc/fastp", paste0(x, ".json")))
    j$summary$before_filtering$total_reads - j$summary$after_filtering$total_reads })
  unal <- sapply(budget$sample, function(x) {
    t <- readLines(file.path(root, "qc", paste0("BOM_ss5__", x, ".hisat2.txt")))
    as.numeric(sub(".*Aligned 0 time: *([0-9]+).*", "\\1", grep("Aligned 0 time", t, value = TRUE)[1])) })

  d <- budget %>%
    mutate(`trimmed`        = as.numeric(trim[sample]),
           `unaligned`      = as.numeric(unal[sample]),
           `UMI duplicate`  = dup$primary[match(sample, dup$sample)] - dup$dedup[match(sample, dup$sample)],
           `multi-mapped`   = get("Unassigned_MultiMapping")[sample],
           `no feature`     = get("Unassigned_NoFeatures")[sample] + get("Unassigned_Ambiguity")[sample],
           # featureCounts Assigned = rRNA + tRNA + mRNA; the read budget tracks tRNA
           # separately, so recover it here or the fate will not reconcile with raw
           tRNA             = get("Assigned")[sample] - rRNA - mRNA) %>%
    select(well, tissue, raw, rRNA, tRNA, mRNA, `multi-mapped`, `no feature`,
           `UMI duplicate`, unaligned, trimmed)
  stopifnot(max(abs(rowSums(d[, 4:11]) - d$raw)) < 5)   # fate must reconcile with raw

  long <- d %>% select(-raw) %>% pivot_longer(-c(well, tissue), names_to = "fate", values_to = "n") %>%
    group_by(well) %>% mutate(frac = n / sum(n)) %>% ungroup()
  lv <- c("rRNA", "tRNA", "mRNA", "multi-mapped", "no feature", "UMI duplicate", "unaligned", "trimmed")
  long$fate <- factor(long$fate, levels = lv)
  pal <- c("rRNA" = "#B2182B", "tRNA" = "#E08214", "mRNA" = "#2166AC", "multi-mapped" = "#F4A582",
           "no feature" = "#D1E5F0", "UMI duplicate" = "#9E9AC8",
           "unaligned" = "#BDBDBD", "trimmed" = "#EEEEEE")

  a <- ggplot(long, aes(well, frac, fill = fate)) +
    geom_col(width = 0.82) +
    facet_grid(~tissue, scales = "free_x", space = "free_x",
               labeller = labeller(tissue = function(x) gsub(" mycelium", "\nmycelium", x))) +
    scale_fill_manual(values = pal, name = NULL) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
    labs(x = NULL, y = "Fraction of raw reads", tag = "a") +
    theme_npj() + theme(legend.position = "right",
                        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5))

  # directly measured rRNA content: subsampled reads mapped against the rRNA loci themselves,
  # independent of the annotation used for counting
  rr <- sapply(budget$sample, function(x) {
    t <- readLines(file.path(root, "qc/rrna", paste0(x, ".summary.txt")))
    g <- function(rx) as.numeric(sub(rx, "\\1", grep(rx, t, value = TRUE)[1]))
    100 * (g(".*Total reads: *([0-9]+).*") - g(".*Aligned 0 time: *([0-9]+).*")) / g(".*Total reads: *([0-9]+).*") })
  bb <- budget %>% mutate(rrna_pct = as.numeric(rr[sample]))
  b <- ggplot(bb, aes(tissue, rrna_pct, colour = tissue)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.3) +
    geom_jitter(width = 0.12, size = 1.1) +
    scale_colour_manual(values = TIS, guide = "none") +
    scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
    ylim(0, 100) +
    labs(x = NULL, y = "Reads matching rRNA loci (%)", tag = "b") + theme_npj()

  c_ <- ggplot(budget, aes(well, mRNA, fill = tissue)) +
    geom_col(width = 0.8) +
    geom_hline(yintercept = 15000, linetype = "22", linewidth = 0.35, colour = "grey20") +
    annotate("text", x = 0.6, y = 15000 * 1.6, label = "retention threshold",
             size = 1.9, hjust = 0, fontface = "italic") +
    scale_fill_manual(values = TIS, name = NULL) +
    scale_y_log10(labels = label_number(big.mark = ","),
                  limits = c(1000, NA), oob = scales::squish) +
    labs(x = NULL, y = "Assigned mRNA counts", tag = "c") +
    theme_npj() + theme(legend.position = "bottom",
                        legend.text = element_text(size = 5),
                        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5))

  save_fig(a / (b | c_) + plot_layout(heights = c(1.1, 1)), "fig1_read_fate", W2, MM(118))
}

# ------------------------------------------------------------------ Figure 2
fig2 <- function() {
  f <- list.files(file.path(root, "qc/testmap"), "\\.summary\\.txt$", full.names = TRUE)
  parse1 <- function(p) {
    t <- readLines(p); g <- function(rx) as.numeric(sub(rx, "\\1", grep(rx, t, value = TRUE)[1]))
    tot <- g(".*Total reads: *([0-9]+).*")
    data.frame(ref = sub("__.*", "", basename(p)),
               unique = 100 * g(".*Aligned 1 time: *([0-9]+).*") / tot,
               multi  = 100 * g(".*Aligned >1 times: *([0-9]+).*") / tot,
               unaligned = 100 * g(".*Aligned 0 time: *([0-9]+).*") / tot)
  }
  tm <- bind_rows(lapply(f, parse1)) %>% group_by(ref) %>%
    summarise(across(everything(), mean), .groups = "drop") %>%
    pivot_longer(-ref, names_to = "cat", values_to = "pct")
  tm$ref <- factor(tm$ref, levels = c("PC9", "PC9.15", "BOM_ss14", "BOM_ss5"))
  tm$cat <- factor(tm$cat, levels = c("unique", "multi", "unaligned"))

  a <- ggplot(tm, aes(ref, pct, fill = cat)) +
    geom_col(width = 0.7) +
    geom_text(data = subset(tm, cat == "unique"),
              aes(label = sprintf("%.1f%%", pct)), vjust = -0.5, size = 2, colour = "#2166AC") +
    scale_fill_manual(values = c(unique = "#2166AC", multi = "#F4A582", unaligned = "#BDBDBD"),
                      name = NULL) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = "% of reads", tag = "a",
         title = "Reference candidates") + theme_npj() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))

  bo <- read_csv(file.path(root, "results/read_budget_BOM_ss5.csv")) %>% select(sample, well, tissue, BOM_ss5 = mRNA)
  pc <- read_csv(file.path(root, "results/read_budget_PC9.15.csv")) %>% select(sample, PC9.15 = mRNA)
  # Plot the per-library PERCENT gain, not the paired absolute counts: the gain is ~7%, which
  # on a log count axis renders as two flat parallel lines and hides the actual result.
  pr <- left_join(bo, pc, by = "sample") %>%
    mutate(gain = 100 * (BOM_ss5 - PC9.15) / PC9.15,
           tissue = factor(tissue, levels = names(TIS))) %>%
    arrange(gain) %>% mutate(well = factor(well, levels = well))
  b <- ggplot(pr, aes(gain, well, colour = tissue)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_segment(aes(x = 0, xend = gain, yend = well), linewidth = 0.3) +
    geom_point(size = 1.3) +
    annotate("text", x = Inf, y = 1.2, hjust = 1.05, size = 2.1, fontface = "italic",
             label = sprintf("mean +%.1f%%", mean(pr$gain))) +
    scale_colour_manual(values = TIS, name = NULL) +
    scale_x_continuous(labels = function(x) paste0("+", x, "%")) +
    labs(x = "Gain in assigned counts, BOM_ss5 vs PC9.15", y = NULL, tag = "b",
         title = "All 16 libraries gain") + theme_npj() +
    theme(legend.position = "bottom", legend.text = element_text(size = 5),
          axis.text.y = element_text(size = 5))

  cov <- read_csv(file.path(root, "results/figure_data/mito_coverage.csv"))
  feat <- read_csv(file.path(root, "results/figure_data/mito_features.csv"))
  blocks <- read.delim(file.path(root, "refs/rRNA_regions_BOM_ss5.tsv"), comment.char = "#", header = FALSE,
                       col.names = c("contig", "start", "end", "name")) %>% filter(contig == "CM148777.1")
  cc <- ggplot(cov, aes(pos / 1000, depth + 1)) +
    geom_rect(data = blocks, inherit.aes = FALSE,
              aes(xmin = start / 1000, xmax = end / 1000, ymin = 1, ymax = Inf),
              fill = "#B2182B", alpha = 0.13) +
    geom_area(fill = "grey25", linewidth = 0) +
    geom_segment(data = feat, inherit.aes = FALSE,
                 aes(x = start / 1000, xend = end / 1000, y = 0.45, yend = 0.45),
                 linewidth = 1.6, colour = "#2166AC") +
    scale_y_log10(labels = label_number(big.mark = ","), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0.35, NA)) +
    labs(x = "Position on mitochondrion CM148777.1 (kb)", y = "Depth + 1", tag = "c",
         title = "Read pile-up falls outside annotated features (blue)") +
    theme_npj()

  save_fig((a | b) / cc + plot_layout(heights = c(1, 0.95)), "fig2_reference_rdna", W2, MM(105))
}

# ------------------------------------------------------------------ Figure 3
fig3 <- function() {
  rd <- function(p, lab) read_csv(p) %>% mutate(set = lab)
  all16 <- rd(file.path(root, "results/dge_BOM_ss5/PCA_rRNArm.csv"), "All 16 libraries")
  keep12 <- rd(file.path(root, "results/dge_BOM_ss5_filtered/PCA_rRNArm.csv"), "12 retained")
  d <- bind_rows(all16, keep12)
  d$Tissue <- gsub("\\.", " ", d$Tissue)
  d$Tissue <- factor(d$Tissue, levels = names(TIS))
  d$set <- factor(d$set, levels = c("All 16 libraries", "12 retained"))
  d$well <- sub(".*_sample_", "", d$sample)
  cors <- d %>% group_by(set) %>% summarise(r = cor(PC1, log10(mrna + 1)), .groups = "drop")

  ab <- ggplot(d, aes(PC1, PC2, colour = Tissue)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey85") +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey85") +
    geom_point(aes(size = log10(mrna + 1))) +
    geom_text(aes(label = well), size = 1.7, vjust = -1.1, show.legend = FALSE) +
    geom_text(data = cors, inherit.aes = FALSE, aes(x = -Inf, y = -Inf,
              label = sprintf("cor(PC1, depth) = %+.3f", r)),
              hjust = -0.06, vjust = -0.8, size = 2.1, fontface = "bold") +
    facet_wrap(~set, scales = "free") +
    scale_colour_manual(values = TIS, name = NULL) +
    scale_size_continuous(range = c(0.8, 3), guide = "none") +
    labs(tag = "a") + theme_npj() + theme(legend.position = "bottom")

  gains <- data.frame(
    metric = factor(rep(c("Genes passing filter", "Genes with tissue effect (LRT)"), each = 2),
                    levels = c("Genes passing filter", "Genes with tissue effect (LRT)")),
    set = factor(rep(c("All 16 libraries", "12 retained"), 2), levels = levels(d$set)),
    n = c(1199, 1666, 82, 281))
  cc <- ggplot(gains, aes(set, n, fill = set)) +
    geom_col(width = 0.6) + geom_text(aes(label = n), vjust = -0.35, size = 2.1) +
    facet_wrap(~metric, scales = "free_y") +
    scale_fill_manual(values = c("All 16 libraries" = "#BDBDBD", "12 retained" = "#2166AC"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = "Genes", tag = "b") + theme_npj()

  save_fig(ab / cc + plot_layout(heights = c(1.45, 1)), "fig3_pca_retention", W2, MM(105))
}

# ------------------------------------------------------------------ Figure 4
fig4 <- function() {
  cr <- as.matrix(read.csv(file.path(root, "results/wgcna/module_trait_correlation.csv"), row.names = 1, check.names = FALSE))
  pv <- as.matrix(read.csv(file.path(root, "results/wgcna/module_trait_pvalue.csv"), row.names = 1, check.names = FALSE))
  fdr <- matrix(p.adjust(pv, "BH"), nrow = nrow(pv), dimnames = dimnames(pv))
  mods <- read.csv(file.path(root, "results/wgcna/gene_modules.csv"))
  sizes <- as.data.frame(table(mods$module), stringsAsFactors = FALSE); names(sizes) <- c("module", "n")

  long <- expand.grid(module = rownames(cr), tissue = colnames(cr), stringsAsFactors = FALSE)
  long$r <- as.vector(cr); long$fdr <- as.vector(fdr)
  long$module <- sub("^ME", "", long$module)
  long <- left_join(long, sizes, by = "module") %>% filter(!is.na(n), module != "grey")
  keep <- long %>% group_by(module) %>% summarise(best = min(fdr), .groups = "drop") %>%
    arrange(best) %>% slice_head(n = 16) %>% pull(module)
  long <- long %>% filter(module %in% keep)
  long$module <- factor(long$module, levels = rev(keep))
  long$star <- ifelse(long$fdr < 0.05, "*", "")

  a <- ggplot(long, aes(tissue, module, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = star), size = 3, vjust = 0.78) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         limits = c(-1, 1), name = "r") +
    scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
    labs(x = NULL, y = NULL, tag = "a",
         title = "Module-tissue correlation (* FDR < 0.05)") + theme_npj()

  conf <- read_csv(file.path(root, "results/wgcna/rrna_confound_test.csv")) %>%
    filter(fdr_tissue < 0.05, abs(r_tissue) > 0.7)
  # keep the tissue label distinct from the correlation columns: selecting `tissue` and
  # renaming r_tissue to `tissue` in the same call produces duplicate names
  cf <- conf %>%
    transmute(module, tis = tissue,
              `with tissue` = r_tissue,
              `with rRNA fraction` = r_rRNAfrac,
              `tissue | rRNA` = r_tissue_partial) %>%
    pivot_longer(c(`with tissue`, `with rRNA fraction`, `tissue | rRNA`),
                 names_to = "what", values_to = "r")
  cf$what <- factor(cf$what, levels = c("with tissue", "with rRNA fraction", "tissue | rRNA"))
  cf$lab <- paste0(cf$module, "\n(", cf$tis, ")")
  b <- ggplot(cf, aes(r, lab, fill = what)) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey70") +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    scale_fill_manual(values = c("with tissue" = "#2166AC", "with rRNA fraction" = "#BDBDBD",
                                 "tissue | rRNA" = "#B2182B"), name = NULL) +
    xlim(-1, 1) + labs(x = "Correlation", y = NULL, tag = "b",
                       title = "Associations survive rRNA control") +
    theme_npj() + theme(legend.position = "bottom", axis.text.y = element_text(size = 5))

  save_fig((a | b) + plot_layout(widths = c(1, 1.15)), "fig4_wgcna", W2, MM(95))
}

cat("writing figures to", FIG, "\n")
fig1(); fig2(); fig3(); fig4()
cat("done\n")

# ============================== SUPPLEMENTARY FIGURES ==============================

figS1 <- function() {   # WGCNA soft-threshold diagnostics
  st <- read_csv(file.path(root, "results/wgcna/soft_threshold.csv"))
  names(st)[1:3] <- c("Power", "SFT.R.sq", "slope")
  a <- ggplot(st, aes(Power, SFT.R.sq)) +
    geom_hline(yintercept = 0.9, linetype = "22", colour = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = 18, colour = "#B2182B", linewidth = 0.4) +
    geom_line(linewidth = 0.3) + geom_point(size = 1) +
    annotate("text", x = 18, y = 0.15, label = "chosen (18)", size = 2, hjust = -0.1, colour = "#B2182B") +
    labs(x = "Soft-threshold power", y = expression("Scale-free topology"~R^2), tag = "a") + theme_npj()
  b <- ggplot(st, aes(Power, mean.k.)) + geom_line(linewidth = 0.3) + geom_point(size = 1) +
    geom_vline(xintercept = 18, colour = "#B2182B", linewidth = 0.4) +
    scale_y_log10() + labs(x = "Soft-threshold power", y = "Mean connectivity", tag = "b") + theme_npj()
  save_fig(a | b, "figS1_soft_threshold", W2, MM(62))
}

figS2 <- function() {   # 3' extension made almost no difference
  rd <- function(p) {
    s <- read.delim(p, check.names = FALSE); rownames(s) <- s$Status
    tot <- colSums(s[, -1]); 100 * as.numeric(s["Assigned", -1]) / tot }
  before <- rd(file.path(root, "qc/fcmp/PC9.15.counts.txt.summary"))
  after  <- rd(file.path(root, "qc/fcmp/PC9.15_3p.counts.txt.summary"))
  d <- data.frame(sample = rep(seq_along(before), 2),
                  pct = c(before, after),
                  set = rep(c("Original annotation", "3' extended (+323 bp median)"), each = length(before)))
  d$set <- factor(d$set, levels = unique(d$set))
  ggp <- ggplot(d, aes(set, pct, group = sample)) +
    geom_line(linewidth = 0.3, colour = "grey60") + geom_point(size = 1.4, colour = "#2166AC") +
    annotate("text", x = 1.5, y = max(d$pct) * 1.06, size = 2.2, fontface = "italic",
             label = sprintf("mean %.1f%% -> %.1f%%", mean(before), mean(after))) +
    labs(x = NULL, y = "Reads assigned to genes (%)",
         title = "Extending 3' ends does not recover the missing reads") + theme_npj()
  save_fig(ggp, "figS2_utr_extension", W1, MM(62))
}

figS3 <- function() {   # tau distribution and marker robustness
  tau <- read_csv(file.path(root, "results/tissue_models_BOM_ss5/tau_specificity.csv"))
  a <- ggplot(tau, aes(tau)) +
    geom_histogram(bins = 40, fill = "#2166AC", colour = "white", linewidth = 0.15) +
    geom_vline(xintercept = 0.85, colour = "#B2182B", linetype = "22", linewidth = 0.4) +
    annotate("text", x = 0.85, y = Inf, label = " tissue-specific", size = 2,
             hjust = 0, vjust = 1.6, colour = "#B2182B") +
    labs(x = expression(tau~"specificity index"), y = "Genes", tag = "a") + theme_npj()
  rb <- data.frame(tissue = c("Exudophore", "Nodule", "Exuding mycelium", "Fuzzy mycelium"),
                   n = c(2, 3, 3, 4), supported = c(37, 22, 1, 1))
  rb$tissue <- factor(rb$tissue, levels = names(TIS))
  b <- ggplot(rb, aes(tissue, supported, fill = tissue)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = paste0(supported, "/50\n(n=", n, ")")), vjust = -0.25, size = 2) +
    scale_fill_manual(values = TIS, guide = "none") +
    scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
    scale_y_continuous(limits = c(0, 50), expand = expansion(mult = c(0, 0.22))) +
    labs(x = NULL, y = "Markers supported by all replicates", tag = "b") + theme_npj()
  save_fig(a | b, "figS3_tau_robustness", W2, MM(64))
}

figS4 <- function() {   # CAZy classes and secretome
  cz <- read.delim(file.path(root, "results/annotation/bom_ss5_cazymes.tsv"))
  d <- as.data.frame(table(cz$cazy_class)); names(d) <- c("class", "n")
  d$class <- factor(d$class, levels = d$class[order(-d$n)])
  a <- ggplot(d, aes(class, n, fill = class)) + geom_col(width = 0.7) +
    geom_text(aes(label = n), vjust = -0.3, size = 2.1) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "CAZy class (Pfam-derived)", y = "Proteins", tag = "a") + theme_npj()
  sec <- read.delim(file.path(root, "results/annotation/bom_ss5_secretome.tsv"))
  tot <- 12521
  d2 <- data.frame(cat = c("Secreted (signal peptide,\nno TM)", "Other annotated", "No Swiss-Prot hit"),
                   n = c(sum(sec$secreted == 1), nrow(sec) - sum(sec$secreted == 1), tot - nrow(sec)))
  d2$cat <- factor(d2$cat, levels = d2$cat)
  b <- ggplot(d2, aes(x = "", y = n, fill = cat)) + geom_col(width = 0.65) +
    geom_text(aes(label = comma(n)), position = position_stack(vjust = 0.5), size = 2.1) +
    scale_fill_manual(values = c("#D55E00", "#2166AC", "#BDBDBD"), name = NULL) +
    coord_flip() + labs(x = NULL, y = "Proteins", tag = "b",
                        title = paste0("Proteome (n = ", comma(tot), ")")) +
    theme_npj() + theme(legend.position = "bottom", legend.text = element_text(size = 5))
  save_fig(a / b + plot_layout(heights = c(1.5, 1)), "figS4_cazy_secretome", W1 * 1.5, MM(85))
}

figS5 <- function() {   # GEM reconstruction stages
  g <- read_csv(file.path(root, "results/figure_data/gem_stats.csv"))
  g$stage <- factor(g$stage, levels = c("draft", "medium", "gapfilled"),
                    labels = c("EC-mapped draft", "+ transport,\nmedium-constrained", "+ targeted\ngapfilling"))
  d <- g %>% select(stage, carrying, blocked) %>%
    pivot_longer(-stage, names_to = "state", values_to = "n")
  d$state <- factor(d$state, levels = c("blocked", "carrying"),
                    labels = c("blocked", "can carry flux"))
  ggp <- ggplot(d, aes(stage, n, fill = state)) +
    geom_col(width = 0.65) +
    geom_text(data = g, inherit.aes = FALSE, aes(stage, carrying, label = comma(carrying)),
              vjust = -0.4, size = 2.1, colour = "#2166AC") +
    scale_fill_manual(values = c("blocked" = "#DDDDDD", "can carry flux" = "#2166AC"), name = NULL) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "Reactions",
         title = "Reactions able to carry flux: 288 -> 2,287") +
    theme_npj() + theme(legend.position = "bottom")
  save_fig(ggp, "figS5_gem_stages", W1 * 1.3, MM(70))
}

figS6 <- function() {   # UMI duplication
  d <- read_csv(file.path(root, "results/figure_data/dedup_rates.csv")) %>%
    mutate(dup = 100 * (primary - dedup) / primary,
           well = sub(".*_sample_", "", sample)) %>%
    left_join(budget %>% select(sample, tissue), by = "sample")
  d$tissue <- factor(d$tissue, levels = names(TIS))
  d <- d %>% arrange(tissue, well); d$well <- factor(d$well, levels = d$well)
  ggp <- ggplot(d, aes(well, dup, fill = tissue)) + geom_col(width = 0.75) +
    scale_fill_manual(values = TIS, name = NULL) +
    labs(x = NULL, y = "UMI duplicate rate (%)",
         title = "PCR/optical duplication per library") +
    theme_npj() + theme(legend.position = "bottom", legend.text = element_text(size = 5),
                        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5))
  save_fig(ggp, "figS6_duplication", W1 * 1.5, MM(62))
}

figS7 <- function() {   # intergenic peaks that located the rDNA
  pk <- read_csv(file.path(root, "results/figure_data/intergenic_peaks.csv")) %>%
    arrange(desc(pct_of_aligned)) %>% slice_head(n = 12) %>%
    mutate(lab = sprintf("%s:%.0f-%.0f kb", sub("\\..*", "", contig), start / 1000, end / 1000),
           rdna = pct_of_aligned > 1)
  pk$lab <- factor(pk$lab, levels = rev(pk$lab))
  ggp <- ggplot(pk, aes(pct_of_aligned, lab, fill = rdna)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", pct_of_aligned)), hjust = -0.15, size = 2) +
    scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#BDBDBD"), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(x = "% of pooled aligned reads", y = NULL,
         title = "Unannotated loci, ranked: rDNA blocks dominate") +
    theme_npj() + theme(axis.text.y = element_text(size = 5))
  save_fig(ggp, "figS7_intergenic_peaks", W1 * 1.6, MM(72))
}

figS1(); figS2(); figS3(); figS4(); figS5(); figS6(); figS7()
cat("supplementary done\n")

# ---------------- Figure for the biology manuscript: exudophore signature ----------------
figExu <- function() {
  norm <- read.csv(file.path(root, "results/dge_BOM_ss5_filtered/normalized_counts_rRNArm.csv"),
                   row.names = 1, check.names = FALSE)
  rb <- read_csv(file.path(root, "results/tissue_models_BOM_ss5/markers_robust.csv")) %>%
    filter(tissue == "Exudophore")
  # markers_robust.csv already carries protein_name; joining the per-tissue marker table
  # duplicated the column as protein_name.x/.y
  d <- rb %>% filter(!is.na(protein_name), protein_name != "") %>%
    arrange(desc(tau)) %>% slice_head(n = 16)

  rs <- read.csv(file.path(root, "metadata/runsheet.csv"), check.names = FALSE)
  tis <- setNames(rs[["Factor Value[Tissue]"]], rs$sample_name)
  sub <- norm[d$gene, , drop = FALSE]
  z <- t(scale(t(as.matrix(sub))))                    # row z-score across libraries
  long <- as.data.frame(z) %>% mutate(gene = rownames(z)) %>%
    pivot_longer(-gene, names_to = "sample", values_to = "z") %>%
    mutate(tissue = factor(tis[sample], levels = names(TIS)),
           well = sub(".*_sample_", "", sample))
  lab <- setNames(paste0(substr(d$protein_name, 1, 42), "  (tau=", sprintf("%.2f", d$tau), ")"), d$gene)
  long$lab <- factor(lab[long$gene], levels = rev(lab[d$gene]))
  long <- long %>% arrange(tissue, well)
  long$well <- factor(long$well, levels = unique(long$well))

  ggp <- ggplot(long, aes(well, lab, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    facet_grid(~tissue, scales = "free_x", space = "free_x",
               labeller = labeller(tissue = function(x) gsub(" ", "\n", x))) +
    scale_fill_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B",
                         midpoint = 0, name = "z-score") +
    labs(x = NULL, y = NULL,
         title = "Replicate-supported exudophore markers") +
    theme_npj() + theme(axis.text.y = element_text(size = 4.6),
                        axis.text.x = element_text(size = 5),
                        legend.position = "right")
  save_fig(ggp, "fig5_exudophore_markers", W2, MM(90))
}
figExu()
