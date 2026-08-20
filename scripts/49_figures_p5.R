#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(circlize); library(ggplot2); library(dplyr); library(tidyr)
  library(readr); library(patchwork); library(scales)
})
options(readr.show_col_types=FALSE)
root <- normalizePath(file.path(dirname(sub("--file=","",grep("--file=",commandArgs(FALSE),value=TRUE)[1])),".."))
FIG <- file.path(root,"latex","figures"); MM <- function(x) x/25.4; W2 <- MM(180)
TIS <- c("Exuding mycelium"="#0072B2","Fuzzy mycelium"="#009E73","Exudophore"="#D55E00","Nodule"="#CC79A7")
FEAT <- c(CDS="#2166AC", rRNA="#B2182B", tRNA="#009E73", intron="#E08214")
theme_npj <- function(base=7) theme_bw(base_size=base) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linewidth=.2,colour="grey92"),
        panel.border=element_rect(linewidth=.3), axis.ticks=element_line(linewidth=.25),
        strip.background=element_rect(fill="grey95",linewidth=.3),
        strip.text=element_text(size=base-.5,face="bold"), legend.key.size=unit(3,"mm"),
        plot.title=element_text(size=base+.5,face="bold"), plot.tag=element_text(size=base+2,face="bold"))
sf <- function(p,n,w,h){ ggsave(file.path(FIG,paste0(n,".pdf")),p,width=w,height=h,units="in",
  device=cairo_pdf,limitsize=FALSE); cat(sprintf("  %-32s %.0f x %.0f mm\n",paste0(n,".pdf"),w*25.4,h*25.4)) }

gff <- read.delim(file.path(root,"refs/BOM_ss5/mitogenome.gff"), comment.char="#", header=FALSE,
                  col.names=c("c","src","type","start","end","score","strand","phase","attr"))
LEN <- 71949

# ---------- F11: circular map (circlize) ----------
f11 <- function() {
  cov <- read_csv(file.path(root,"results/figure_data/mito_coverage.csv"))
  bin <- cov %>% mutate(b=(pos %/% 250)*250) %>% group_by(b) %>%
    summarise(d=mean(depth)+1, .groups="drop")
  dup <- data.frame(s=c(63979,66816,34388,68781), e=c(64450,67285,34717,69108))
  pdf(file.path(FIG,"fig11_circular_map.pdf"), width=MM(120), height=MM(120))
  circos.clear()
  circos.par(start.degree=90, gap.degree=4, cell.padding=c(0,0,0,0), track.margin=c(0.006,0.006))
  circos.initialize(factors="mt", xlim=c(0,LEN))
  # coverage
  circos.track(factors="mt", ylim=c(0,log10(max(bin$d))), track.height=0.20, bg.border=NA,
    panel.fun=function(x,y) circos.lines(bin$b, log10(bin$d), type="h", col="grey35", lwd=0.35))
  circos.yaxis(side="left", at=c(0,2,4,6), labels=c("1","100","10k","1M"),
               labels.cex=0.3, sector.index="mt")
  # features, forward and reverse
  for (st in c("+","-")) {
    circos.track(factors="mt", ylim=c(0,1), track.height=0.075, bg.border=NA, panel.fun=function(x,y){
      d <- gff[gff$strand==st & gff$type %in% names(FEAT),]
      if (nrow(d)) circos.rect(d$start, 0.15, d$end, 0.85, col=FEAT[d$type], border=NA)
    })
  }
  # introns (strandless)
  circos.track(factors="mt", ylim=c(0,1), track.height=0.06, bg.border=NA, panel.fun=function(x,y){
    d <- gff[gff$type=="intron",]
    if (nrow(d)) circos.rect(d$start,0.2,d$end,0.8, col=FEAT["intron"], border=NA)
  })
  # axis
  circos.track(factors="mt", ylim=c(0,1), track.height=0.05, bg.border=NA, panel.fun=function(x,y)
    circos.axis(h="top", major.at=seq(0,70000,10000), labels=paste0(seq(0,70,10),"kb"),
                labels.cex=0.35, major.tick.length=0.3))
  # duplication links
  for (i in seq(1,nrow(dup),2))
    circos.link("mt", c(dup$s[i],dup$e[i]), "mt", c(dup$s[i+1],dup$e[i+1]),
                col="#88000055", border=NA)
  text(0,0.10,"Pleurotus ostreatus", cex=0.42, font=3)
  text(0,0.02,"BOM_ss5 mitochondrion", cex=0.42)
  text(0,-0.06,"71,949 bp", cex=0.38, col="grey30")
  legend("bottomleft", legend=c(names(FEAT),"duplication"), fill=c(FEAT,"#88000055"),
         border=NA, bty="n", cex=0.34)
  dev.off(); circos.clear()
  cat("  fig11_circular_map.pdf            120 x 120 mm\n")
}

# ---------- F12: comparative size vs intron content ----------
f12 <- function() {
  inv <- read_csv(file.path(root,"results/mito_comparative/inventory.csv"))
  intr <- lapply(inv$acc, function(a){
    p <- file.path(root,"results/mito_comparative",paste0(a,".introns.tbl"))
    if (!file.exists(p)) return(data.frame(acc=a,n=0,span=0))
    x <- readLines(p); x <- x[!startsWith(x,"#")]
    if (!length(x)) return(data.frame(acc=a,n=0,span=0))
    f <- strsplit(trimws(x),"\\s+")
    s <- sapply(f,function(z) abs(as.numeric(z[9])-as.numeric(z[8]))+1)
    data.frame(acc=a,n=length(s),span=sum(s))
  }) %>% bind_rows()
  d <- inv %>% left_join(intr, by="acc") %>%
    mutate(lab=sub("\\..*","",acc),
           grp=ifelse(acc %in% c("CM148777.1","CM148778.1"),"this study",
               ifelse(grepl("^CM057219|^OR0301|^PX7243", acc),"P. ostreatus (other)","other Pleurotus")))
  r <- cor(d$length, d$span)
  a <- ggplot(d, aes(span, length/1000, colour=grp)) +
    geom_smooth(method="lm", se=FALSE, colour="grey60", linewidth=.35, formula=y~x) +
    geom_point(size=1.6) +
    ggrepel::geom_text_repel(aes(label=lab), size=1.8, max.overlaps=20, seed=1) +
    scale_colour_manual(values=c("this study"="#D55E00","P. ostreatus (other)"="#2166AC",
                                 "other Pleurotus"="grey55"), name=NULL) +
    labs(x="Group I intron span detected (bp)", y="Mitogenome size (kb)", tag="a",
         title=sprintf("Intron content tracks genome size (r = %+.2f)", r)) +
    theme_npj() + theme(legend.position="bottom")
  b <- ggplot(d, aes(orfs_ge100aa, length/1000, colour=grp)) +
    geom_smooth(method="lm", se=FALSE, colour="grey60", linewidth=.35, formula=y~x) +
    geom_point(size=1.6) +
    scale_colour_manual(values=c("this study"="#D55E00","P. ostreatus (other)"="#2166AC",
                                 "other Pleurotus"="grey55"), guide="none") +
    labs(x="ORFs (>=100 aa)", y=NULL, tag="b",
         title=sprintf("and ORF count more so (r = %+.2f)", cor(d$length,d$orfs_ge100aa))) +
    theme_npj()
  sf(a|b, "fig12_comparative_size", W2, MM(80))
}

# ---------- F13: functional evidence for the uncharacterised ORFs ----------
f13 <- function() {
  cons <- read_csv(file.path(root,"results/mito/orf_conservation.csv"))
  d <- cons %>% mutate(status=case_when(
        genomes>=6 & libs>=8 & in_rrna==0 ~ "conserved + expressed",
        genomes>=6 ~ "conserved, low expression",
        genomes<=1 ~ "strain-specific", TRUE ~ "partial"))
  a <- ggplot(d, aes(genomes, libs, colour=status, size=reads+1)) +
    geom_point(alpha=.85) +
    scale_colour_manual(values=c("conserved + expressed"="#B2182B",
      "conserved, low expression"="#2166AC","partial"="grey60","strain-specific"="#E08214"), name=NULL) +
    scale_size_continuous(range=c(.8,4.2), guide="none") +
    scale_x_continuous(breaks=seq(0,11,2)) +
    labs(x="Present in n of 11 other Pleurotus mitogenomes", y="Detected in n of 16 libraries",
         tag="a", title="Uncharacterised ORFs: conservation against expression") +
    theme_npj() + theme(legend.position="bottom", legend.text=element_text(size=5))
  rg <- read_csv(file.path(root,"results/mitocarta/retrograde_per_library.csv"))
  rg$tissue <- factor(rg$tissue, levels=names(TIS))
  b <- ggplot(rg, aes(tissue, ratio, colour=tissue)) +
    geom_boxplot(outlier.shape=NA, width=.5, linewidth=.3) +
    geom_jitter(width=.12, size=1.4) +
    scale_colour_manual(values=TIS, guide="none") +
    scale_x_discrete(labels=function(x) gsub(" ","\n",x)) +
    labs(x=NULL, y="mtDNA-encoded / nuclear mitochondrial", tag="b",
         title="Mito-nuclear balance by tissue") + theme_npj()
  sf(a|b, "fig13_orfs_retrograde", W2, MM(78))
}
cat("P5 figures:\n"); f11(); f12(); f13()
