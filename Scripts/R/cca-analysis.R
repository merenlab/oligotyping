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
csv_cca <- cca(csv[,-1])

C <- function(){
    plot(csv_cca, type = "n")
    title(title_text, col.main = "gray")
    text(csv_cca, dis = "sites", col="red", cex=0.5, labels=rownames(csv))
    points(csv_cca, display = "sp", col = "gray", bg = "red", cex =0.4)
}

C()
tk_messageBox(message="Press a key")

# PDF
pdf_output <- paste(output_file_prefix,".pdf",sep="")
pdf(pdf_output)
C()
sprintf("CCA result PDF: '%s'", pdf_output)
dev.off()

# PNG
png_output <- paste(output_file_prefix,".png",sep="")
png(png_output, width = 1200, height = 1000, units = "px", pointsize = 12, bg = "transparent", type = c("cairo", "cairo-png", "Xlib", "quartz"))
C()
sprintf("CCA result PNG: '%s'", png_output)
dev.off()
