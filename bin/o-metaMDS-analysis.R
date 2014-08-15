#!/usr/bin/env Rscript

library("vegan")
library(gtools)

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


csv <- read.csv(csv_path, header=TRUE, sep="\t")
rownames(csv) <- csv[,1]

fit <- metaMDS(csv[,-1], distance=distance)

P <- function(){
    plot(fit, cex = 0.7, main = title_text, sub = paste("Distance metric: ",distance,sep=""))
    ordilabel (fit, display = c('sites'))
}


# PDF
pdf_output <- paste(output_file_prefix,".pdf",sep="")
pdf(pdf_output, width = 8, height = 8, pointsize = 8, family='Helvetica')
P()
sprintf("Clustering result PDF: '%s'", pdf_output)
dev.off()

# PNG
png_output <- paste(output_file_prefix,".png",sep="")
png(png_output, width = 800, height = 600, units = "px", pointsize = 12, bg = "transparent", type = c("cairo", "cairo-png", "Xlib", "quartz"))
P()
sprintf("Clustering result PNG: '%s'", png_output)
dev.off()
