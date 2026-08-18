#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(patchwork); library(scales)
})
options(readr.show_col_types = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1])), ".."))
FIG <- file.path(root,"latex","figures"); MM <- function(x) x/25.4; W2 <- MM(180)
theme_npj <- function(base=7) theme_bw(base_size=base) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linewidth=.2,colour="grey92"),
        panel.border=element_rect(linewidth=.3), axis.ticks=element_line(linewidth=.25),
        strip.background=element_rect(fill="grey95",linewidth=.3),
        strip.text=element_text(size=base-.5,face="bold"), legend.key.size=unit(3,"mm"),
        plot.title=element_text(size=base+.5,face="bold"), plot.tag=element_text(size=base+2,face="bold"))
sf <- function(p,n,w,h){ ggsave(file.path(FIG,paste0(n,".pdf")),p,width=w,height=h,units="in",
  device=cairo_pdf,limitsize=FALSE); cat(sprintf("  %-30s %.0f x %.0f mm\n",paste0(n,".pdf"),w*25.4,h*25.4)) }

# ---------- F8: mitogenome map, stock vs corrected, with coverage ----------
f8 <- function() {
  cov <- read_csv(file.path(root,"results/figure_data/mito_coverage.csv"))
  gff <- read.delim(file.path(root,"refs/BOM_ss5/mitogenome.gff"), comment.char="#", header=FALSE,
                    col.names=c("c","src","type","start","end","score","strand","phase","attr"))
  stock <- read.delim(file.path(root,"refs/BOM_ss5/BOM_ss5_genomic.gtf"), header=FALSE, comment.char="#",
                      col.names=c("c","src","type","start","end","score","strand","phase","attr"))
  stock <- stock %>% filter(c=="CM148777.1", type=="gene")
  a <- ggplot(cov, aes(pos/1000, depth+1)) + geom_area(fill="grey30", linewidth=0) +
    scale_y_log10(labels=comma, expand=c(0,0)) +
    labs(x=NULL, y="Depth + 1", tag="a", title="Pooled read depth across the mitochondrion") +
    theme_npj()
  lv <- c("CDS","rRNA","tRNA","intron")
  gg <- gff %>% filter(type %in% lv) %>% mutate(track="Corrected (this study)")
  ss <- stock %>% mutate(type=ifelse(grepl('protein_coding',attr),"CDS","tRNA"), track="Stock annotation")
  d <- bind_rows(gg %>% select(start,end,type,track), ss %>% select(start,end,type,track))
  d$type <- factor(d$type, levels=lv)
  d$track <- factor(d$track, levels=c("Stock annotation","Corrected (this study)"))
  b <- ggplot(d) +
    geom_rect(aes(xmin=start/1000, xmax=end/1000, ymin=0, ymax=1, fill=type)) +
    facet_wrap(~track, ncol=1, strip.position="left") +
    scale_fill_manual(values=c(CDS="#2166AC", rRNA="#B2182B", tRNA="#009E73", intron="#E08214"), name=NULL) +
    scale_y_continuous(breaks=NULL, expand=c(0,0)) +
    labs(x="Position on CM148777.1 (kb)", y=NULL, tag="b") +
    theme_npj() + theme(legend.position="bottom", strip.placement="outside",
                        strip.text.y.left=element_text(angle=0, size=5.5))
  sf(a/b + plot_layout(heights=c(1,0.85)), "fig8_mitogenome_map", W2, MM(78))
}

# ---------- F9: intron catalogue ----------
f9 <- function() {
  ic <- read_csv(file.path(root,"results/mito/intron_catalogue.csv")) %>%
    filter(!is.na(mean_efficiency)) %>%
    mutate(real = mean_efficiency > 0.5,
           rfam_hit = nchar(rfam) > 0,
           in_rrna = start >= 1000 & start <= 7200)
  a <- ggplot(ic, aes(length, mean_efficiency, colour=real, size=total_reads)) +
    geom_hline(yintercept=0.5, linetype="22", linewidth=.3, colour="grey50") +
    geom_point(alpha=.75) +
    scale_colour_manual(values=c(`TRUE`="#B2182B",`FALSE`="grey65"),
                        labels=c("artefact (<50%)","spliced intron"), name=NULL) +
    scale_size_continuous(range=c(.6,4), guide="none") + scale_x_log10(labels=comma) +
    labs(x="Junction length (bp)", y="Mean splicing efficiency", tag="a",
         title="Six junctions splice; fifty do not") + theme_npj() + theme(legend.position="bottom")
  cn <- ic %>% count(in_rrna, real) %>%
    mutate(grp=ifelse(in_rrna,"inside mt-rRNA block","elsewhere"),
           lab=ifelse(real,"spliced intron","artefact"))
  b <- ggplot(cn, aes(grp, n, fill=lab)) + geom_col(width=.6) +
    geom_text(aes(label=n), position=position_stack(vjust=.5), size=2.2, colour="white") +
    scale_fill_manual(values=c(`spliced intron`="#B2182B",`artefact`="grey65"), name=NULL) +
    labs(x=NULL, y="Junctions", tag="b",
         title="Artefacts concentrate in the\nhigh-coverage rRNA region") +
    theme_npj() + theme(legend.position="bottom")
  real <- ic %>% filter(real) %>% arrange(desc(total_reads)) %>%
    mutate(lab=sprintf("%d-%d (%d bp)%s", start, end, length, ifelse(rfam_hit," *","")))
  real$lab <- factor(real$lab, levels=rev(real$lab))
  c_ <- ggplot(real, aes(total_reads, lab)) +
    geom_col(fill="#B2182B", width=.65) +
    geom_text(aes(label=sprintf("%d/16 libs", libraries)), hjust=-0.1, size=2) +
    scale_x_log10(labels=comma, expand=expansion(mult=c(0,.35))) +
    labs(x="Supporting split reads", y=NULL, tag="c",
         title="Confirmed introns (* Rfam group I)") +
    theme_npj() + theme(axis.text.y=element_text(size=5.5))
  sf((a|b)/c_ + plot_layout(heights=c(1,.8)), "fig9_introns", W2, MM(100))
}

# ---------- F10: expressed mitotranscriptome and candidate ORFs ----------
f10 <- function() {
  fc <- read.delim(file.path(root,"results/mito/mito_counts.txt"), comment.char="#", check.names=FALSE)
  gff <- readLines(file.path(root,"refs/BOM_ss5/mitogenome.gff"))
  gff <- gff[!startsWith(gff,"#")]
  ids <- sub(".*ID=([^;]+);.*","\\1",gff); nms <- sub(".*Name=([^;]+);.*","\\1",gff)
  typ <- vapply(strsplit(gff,"\t"), `[`, "", 3)
  map <- setNames(nms, ids); tmap <- setNames(typ, ids)
  cnt <- data.frame(id=fc$Geneid, total=rowSums(fc[,7:ncol(fc)]))
  cnt$name <- map[cnt$id]; cnt$type <- tmap[cnt$id]
  cnt <- cnt %>% filter(total>0, type %in% c("CDS","rRNA")) %>% arrange(desc(total)) %>% head(16)
  cnt$name <- factor(cnt$name, levels=rev(cnt$name))
  a <- ggplot(cnt, aes(total, name, fill=type)) + geom_col(width=.7) +
    scale_fill_manual(values=c(CDS="#2166AC", rRNA="#B2182B"), name=NULL) +
    scale_x_log10(labels=comma, expand=expansion(mult=c(0,.1))) +
    labs(x="Pooled reads", y=NULL, tag="a", title="Expressed mitochondrial features") +
    theme_npj() + theme(legend.position="bottom", axis.text.y=element_text(size=5))
  orf <- read_csv(file.path(root,"results/mito/candidate_orfs.csv")) %>%
    mutate(cls=ifelse(in_rRNA_block==1,"inside rRNA block (discounted)","credible"))
  b <- ggplot(orf, aes(start/1000, total_reads+1, colour=cls, size=libraries_detected)) +
    geom_point(alpha=.8) +
    scale_colour_manual(values=c(credible="#2166AC",`inside rRNA block (discounted)`="grey65"), name=NULL) +
    scale_size_continuous(range=c(.7,3.4), name="libraries") +
    scale_y_log10(labels=comma) +
    labs(x="Position (kb)", y="Reads + 1", tag="b",
         title="Uncharacterised ORFs: expression by position") +
    theme_npj() + theme(legend.position="bottom", legend.box="vertical",
                        legend.text=element_text(size=5))
  sf(a|b, "fig10_mito_expression", W2, MM(80))
}
cat("P4 figures:\n"); f8(); f9(); f10()
