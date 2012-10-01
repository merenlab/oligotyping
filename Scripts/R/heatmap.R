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
    title_text <- "Unknown Title"

if(invalid(output_file_prefix))
    output_file_prefix <- "unknown"

csv <- read.csv(csv_path, header=TRUE, sep="\t")
rownames(csv) <- csv[,1]
m = as.matrix(csv[,-1])

my.colors <- colorRampPalette(c("gray10","gray20","gray30","gray40","gray50","gray60","gray80","gray90","gray100"))

hmap <- function(){
    heatmap.2(t(m),
              dendrogram="column",
              scale="row",
              col = my.colors(512),
              density.info="none",
              trace="none",
              cexRow=0.9,
              cexCol=0.9,
              symm=F,
              symkey=T,
              symbreaks=T, 
              labRow=NA,
              margins = c(5,5),
              key=FALSE,
              colsep=1:nrow(csv),
              rowsep=1:ncol(csv),
              sepcolor='white')

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

