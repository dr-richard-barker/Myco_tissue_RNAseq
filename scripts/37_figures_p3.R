#!/usr/bin/env Rscript
# Figures for the systems-biology manuscript (P3).
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(patchwork); library(scales)
})
options(readr.show_col_types = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
FIG <- file.path(root, "latex", "figures"); MM <- function(x) x/25.4; W1 <- MM(88); W2 <- MM(180)
TIS <- c("Exuding mycelium"="#0072B2","Fuzzy mycelium"="#009E73","Exudophore"="#D55E00","Nodule"="#CC79A7")
theme_npj <- function(base=7) theme_bw(base_size=base) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linewidth=.2,colour="grey92"),
        panel.border=element_rect(linewidth=.3), axis.ticks=element_line(linewidth=.25),
        strip.background=element_rect(fill="grey95",linewidth=.3),
        strip.text=element_text(size=base-.5,face="bold"), legend.key.size=unit(3,"mm"),
        plot.title=element_text(size=base+.5,face="bold"), plot.tag=element_text(size=base+2,face="bold"))
save_fig <- function(p,n,w,h){ ggsave(file.path(FIG,paste0(n,".pdf")),p,width=w,height=h,units="in",
  device=cairo_pdf,limitsize=FALSE); cat(sprintf("  %-32s %.0f x %.0f mm\n",paste0(n,".pdf"),w*25.4,h*25.4)) }

# ---- Fig 6: model overview and the presence/expression caveat ----
fig6 <- function() {
  s <- read_csv(file.path(root,"models/tissue/tissue_model_summary.csv"))
  s$tissue <- factor(s$tissue, levels=names(TIS))
  a <- ggplot(s, aes(tissue, reactions, fill=tissue)) + geom_col(width=.65) +
    geom_text(aes(label=reactions), vjust=-.35, size=2.1) +
    scale_fill_manual(values=TIS, guide="none") +
    scale_x_discrete(labels=function(x) gsub(" ","\n",x)) +
    scale_y_continuous(expand=expansion(mult=c(0,.16))) +
    labs(x=NULL,y="Reactions retained",tag="a",title="Context-specific model size") + theme_npj()

  pa <- read_csv(file.path(root,"results/tissue_metabolism/reaction_presence.csv"))
  cnt <- pa %>% mutate(n=Exuding_mycelium+Fuzzy_mycelium+Exudophore+Nodule) %>% count(n)
  b <- ggplot(cnt, aes(factor(n), nn <- n)) + geom_col(fill="#666666", width=.65) +
    geom_text(aes(label=after_stat(y)), stat="identity", vjust=-.35, size=2.1) +
    scale_y_continuous(expand=expansion(mult=c(0,.16))) +
    labs(x="Number of tissues retaining the reaction", y="Reactions", tag="b",
         title="Reaction sharing across models") + theme_npj()

  # the caveat: reactions unique to the exudophore model, against the expression behind them
  re <- read_csv(file.path(root,"results/tissue_metabolism/reaction_expression.csv"))
  uq <- pa %>% filter(Exudophore==1, Exuding_mycelium+Fuzzy_mycelium+Nodule==0) %>% pull(reaction)
  d <- re %>% mutate(cls=ifelse(reaction %in% uq,"unique to exudophore model","shared"))
  c_ <- ggplot(d, aes(cls, ratio, colour=cls)) +
    geom_hline(yintercept=1, linetype="22", linewidth=.3, colour="grey50") +
    geom_jitter(width=.18, size=.5, alpha=.45) +
    geom_boxplot(outlier.shape=NA, fill=NA, width=.45, linewidth=.35) +
    scale_colour_manual(values=c("shared"="#999999","unique to exudophore model"="#D55E00"), guide="none") +
    scale_y_log10() + scale_x_discrete(labels=function(x) gsub(" ","\n",x)) +
    labs(x=NULL, y="Exudophore / highest other tissue", tag="c",
         title="Model membership does not imply\nelevated expression") + theme_npj()
  save_fig((a|b|c_) + plot_layout(widths=c(1,1.1,.95)), "fig6_tissue_models", W2, MM(62))
}

# ---- Fig 7: reaction-level expression contrast ----
fig7 <- function() {
  re <- read_csv(file.path(root,"results/tissue_metabolism/reaction_expression.csv"))
  re <- re %>% mutate(strong = ratio>=2 & Exudophore>=5)
  a <- ggplot(re, aes(max_other+1, Exudophore+1, colour=strong)) +
    geom_abline(slope=1, intercept=0, linewidth=.3, colour="grey60") +
    geom_abline(slope=1, intercept=log10(2), linetype="22", linewidth=.3, colour="#D55E00") +
    geom_point(size=.6, alpha=.5) +
    scale_x_log10(labels=comma) + scale_y_log10(labels=comma) +
    scale_colour_manual(values=c(`TRUE`="#D55E00",`FALSE`="grey70"), guide="none") +
    labs(x="Highest other tissue (normalised)", y="Exudophore (normalised)", tag="a",
         title=sprintf("%d of %d reactions >=2x elevated", sum(re$strong), nrow(re))) + theme_npj()

  # curated enzymes underpinning the proposed biochemistry
  want <- c("Methanol:oxygen"="Alcohol/methanol oxidase",
            "Acetaldehyde:NAD\\+ oxidoreductase"="Aldehyde dehydrogenase",
            "D-mannitol"="Mannitol dehydrogenase",
            "D-sorbitol 2-dehydrogenase"="Sorbitol dehydrogenase",
            "thiosulfate:cyanide"="Rhodanese (cyanide detox)",
            "linoleate:oxygen"="Linoleate 8R-dioxygenase (oxylipin)",
            "formate  dehydrogenase|formate dehydrogenase"="Formate dehydrogenase",
            "cytochrome-c peroxidase"="Cytochrome-c peroxidase")
  rows <- lapply(names(want), function(p){
    h <- re %>% filter(grepl(p, name, ignore.case=TRUE)) %>% arrange(desc(ratio)) %>% slice_head(n=1)
    if(nrow(h)==0) return(NULL)
    h$label <- want[[p]]; h })
  d <- bind_rows(rows)
  if (nrow(d)) {
    long <- d %>% select(label, Exudophore, Exuding_mycelium, Fuzzy_mycelium, Nodule) %>%
      pivot_longer(-label, names_to="tissue", values_to="v") %>%
      mutate(tissue=factor(gsub("_"," ",tissue), levels=names(TIS)))
    long$label <- factor(long$label, levels=rev(d$label[order(d$ratio)]))
    b <- ggplot(long, aes(v+1, label, fill=tissue)) +
      geom_col(position=position_dodge(width=.8), width=.75) +
      scale_fill_manual(values=TIS, name=NULL) + scale_x_log10(labels=comma) +
      labs(x="Normalised expression + 1", y=NULL, tag="b",
           title="Enzymes of the proposed exudophore biochemistry") +
      theme_npj() + theme(legend.position="bottom", legend.text=element_text(size=5),
                          axis.text.y=element_text(size=5.5))
    save_fig(a/b + plot_layout(heights=c(1,1.25)), "fig7_reaction_expression", W2, MM(115))
  }
}
cat("P3 figures:\n"); fig6(); fig7()
