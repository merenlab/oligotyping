#!/usr/bin/env Rscript

X11(width=12, height=10)
library("vegan")
library(tcltk)
library(gtools)
library(gplots)

args <- commandArgs(trailingOnly = TRUE)
csv_path <- args[1]
output_file_prefix <- args[2]
title_text <- args[3]


if(invalid(title_text))
    title_text <- ""

if(invalid(output_file_prefix))
    output_file_prefix <- "unknown"

csv <- read.csv(csv_path, header=TRUE, sep="\t")
rownames(csv) <- csv[,1]
m = as.matrix(csv[,-1])


hmap <- function(){
    heatmap.2(t(m),
              dendrogram="both",
              scale="row",
              trace="none",
              cexRow=0.4,
              cexCol=1.5,
              distfun=function(m) vegdist(m,method="bray"), #"manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "morisita", "horn", "mountford", "raup" , "binomial" or "chao"
              hclustfun=function(x) hclust(x,method="ward"), # "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
              symm=F,
              symkey=T,
              symbreaks=T, 
              col=colorRampPalette(c("red","green","green4","violet","purple"))(100),
              labRow=NA,
              margins = c(10,10),
              )

    title(title_text, col.main = "gray")
}

hmap()
tk_messageBox(message="Press a key")

# PDF
pdf_output <- paste(output_file_prefix,".pdf",sep="")
pdf(pdf_output)
hmap()
sprintf("Heatmap result PDF: '%s'", pdf_output)
dev.off()

# PNG
png_output <- paste(output_file_prefix,".png",sep="")
png(png_output, width = 1200, height = 1000, units = "px", pointsize = 12, bg = "transparent", type = c("cairo", "cairo-png", "Xlib", "quartz"))
sprintf("Heatmap result PNG: '%s'", png_output)
hmap()
dev.off()

