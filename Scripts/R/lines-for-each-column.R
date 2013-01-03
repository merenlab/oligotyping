#!/usr/bin/env Rscript
#
# generates lines for each column in a given TAB delimited file:
#

library(reshape)
library(gtools)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
csv_path <- args[1]
output_file_prefix <- args[2]
title_text <- args[3]

if(invalid(title_text))
    title_text <- "Unknown Title"

if(invalid(output_file_prefix))
    output_file_prefix <- "unknown"

data <- as.data.frame(read.csv(csv_path, header=TRUE, sep="\t"))
row.names <- data$samples
col.names <- colnames(data)

df <- melt(data ,  id = 'samples', variable_name = 'bins')

P <- function(){
    ggplot(df, aes(x = samples, y = value, colour = bins, group = bins)) + geom_line() + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'top') + labs(x='', y='Number of Reads')
}

pdf_w = 5 + (length(row.names) / 10)
png_w = pdf_w * 100

# PDF
pdf_output <- paste(output_file_prefix,".pdf",sep="")
pdf(pdf_output, width = pdf_w, height = 4)
P()
sprintf("Lines PDF: '%s'", pdf_output)

# PNG
png_output <- paste(output_file_prefix,".png",sep="")
png(png_output, width = png_w, height = 400, units = "px", pointsize = 12, bg = "transparent", type = c("cairo", "cairo-png", "Xlib", "quartz"))
P()
sprintf("Lines PNG: '%s'", png_output)
