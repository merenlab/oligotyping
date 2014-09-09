#!/usr/bin/env Rscript

library("vegan")
library(gtools)

# DATA DIST TITLE OUTPUT_PREFIX

args <- commandArgs(trailingOnly = TRUE)
csv_path <- args[1]
distance <- args[2]
title_text <- args[3]
output_file_prefix <- args[4]

if(invalid(distance))
    distance <- "horn"

if(invalid(title_text))
    title_text <- "Unknown Title"

display <- FALSE
if(invalid(output_file_prefix)){
    output_file_prefix <- "unknown"
    display <- TRUE
}

if (display == TRUE){
    X11(width=12, height=10)
}


csv <- read.csv(csv_path, header=TRUE, sep="\t")
rownames(csv) <- csv[,1]

d <- vegdist(csv[,-1], method=distance) #"manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "morisita", "horn", "mountford", "raup" , "binomial" or "chao"
fit <- hclust(d, method="ward") # "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"

num_samples <- length(row.names)
pdf_w <- num_samples / 4
if(num_samples < 32)
    pdf_w <- 8
png_w = pdf_w * 100

P <- function(){
    plot(fit, labels=rownames(csv), cex = 0.7, main = title_text, sub = paste("Distance metric: ",distance,sep=""), xlab = '') # display dendogram
}

if(display == TRUE){
    P()
    tk_messageBox(message="Press a key")
}

# PDF
pdf_output <- paste(output_file_prefix,".pdf",sep="")
pdf(pdf_output, width = pdf_w, pointsize = 8, family = 'Helvetica')
P()
sprintf("Clustering result PDF: '%s'", pdf_output)
dev.off()

# PNG
png_output <- paste(output_file_prefix,".png",sep="")
png(png_output, width = png_w, height = 1000, units = "px", pointsize = 12, bg = "transparent", type = c("cairo", "cairo-png", "Xlib", "quartz"))
P()
sprintf("Clustering result PNG: '%s'", png_output)
dev.off()
