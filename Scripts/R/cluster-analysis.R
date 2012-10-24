#!/usr/bin/env Rscript

X11(width=12, height=10)
library("vegan")
library(tcltk)
library(gtools)

args <- commandArgs(trailingOnly = TRUE)
csv_path <- args[1]
output_file_prefix <- args[2]
title_text <- args[3]


if(invalid(title_text))
    title_text <- "Unknown Title"

if(invalid(output_file_prefix))
    output_file_prefix <- "unknown"

csv <- read.csv(csv_path, header=TRUE, sep="\t")
rownames(csv) <- csv[,1]


d <- vegdist(csv[,-1], method="horn") #"manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "morisita", "horn", "mountford", "raup" , "binomial" or "chao"
fit <- hclust(d, method="ward") # "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"

P <- function(){
    plot(fit, labels=rownames(csv), cex = 0.7) # display dendogram
}

P()
tk_messageBox(message="Press a key")

# PDF
pdf_output <- paste(output_file_prefix,".pdf",sep="")
pdf(pdf_output)
P()
sprintf("Clustering result PDF: '%s'", pdf_output)
dev.off()

# PNG
png_output <- paste(output_file_prefix,".png",sep="")
png(png_output, width = 1200, height = 1000, units = "px", pointsize = 12, bg = "transparent", type = c("cairo", "cairo-png", "Xlib", "quartz"))
P()
sprintf("Clustering result PNG: '%s'", png_output)
dev.off()
