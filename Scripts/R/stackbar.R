#!/usr/bin/env Rscript
#
# generates stack-bar chart from environment file.
#

library(reshape)
library(gtools)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
environment_file_path <- args[1]
title_text <- args[2]
output_file_prefix <- args[3]
colors_file <- args[4]


if(invalid(title_text))
    title_text <- "Unknown Title"

if(invalid(output_file_prefix))
    output_file_prefix <- "unknown"

df <- as.data.frame(read.csv(environment_file_path, header=FALSE, sep="\t"))
names(df) <- c('oligo', 'sample', 'abundance')

df$oligo <- reorder(df$oligo, df$abundance, FUN=sum)
df$oligo <- factor(df$oligo, levels=levels(df$oligo))
df <- df[order(df$oligo), ]

if(!invalid(colors_file)){
    colors_list <- as.data.frame(read.csv(colors_file, header=FALSE, sep="\t"))
    names(colors_list) <- c('oligo', 'color')
    colors_list$oligo <- factor(colors_list$oligo, levels=rev(levels(df$oligo)))
    colors_list <- colors_list[order(colors_list$oligo), ]
}

P <- function(df){
    p <- ggplot(df, aes(x=factor(sample), y=abundance, fill = oligo))
    p <- p + geom_bar(position="fill", stat = "identity", width=0.90, colour = 'black')
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10))
    p <- p + theme(axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(legend.text = element_text(size = 8, family='mono'), legend.position = 'none')
    p <- p + labs(x='', y='Oligotypes', title=paste('Distribution of Oligotypes Among Samples for', title_text, sep=' '))
    p <- p + scale_y_continuous(breaks = NULL)
    p <- p + guides(fill = guide_legend(nrow = 25))
    p <- p + theme(plot.title = element_text(hjust=0, vjust=1))
    p <- p + coord_cartesian(ylim=c(-0.01, 1.01))

    # colors file was provided?
    if(!invalid(colors_file)){
        p <- p + scale_fill_manual(limits = colors_list[[1]], values = as.vector(colors_list[[2]]))
    }

    print(p)
}

num_samples <- length(levels(factor(df$sample)))
pdf_w <- num_samples / 4
if(num_samples < 32)
    pdf_w <- 8

png_w <- pdf_w * 100

# PDF
pdf_output <- paste(output_file_prefix,".pdf",sep="")
pdf(pdf_output, width = pdf_w, height = 8, pointsize = 8, family='Helvetica')
P(df)
sprintf("Stackbar PDF: '%s'", pdf_output)

# PNG
png_output <- paste(output_file_prefix,".png",sep="")
png(png_output, width = png_w, height = 800, units = "px", pointsize = 12, bg = "transparent", type = c("cairo", "cairo-png", "Xlib", "quartz"))
P(df)
sprintf("Stackbar PNG: '%s'", png_output)
