#!/usr/bin/env Rscript
#
# generates stack-bar chart from environment file.
#

suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(reshape))

# command line options
option_list <- list(
		make_option(c("-o", "--output_file_prefix"), default="unknown",
				help = "Output file prefix for visualization files [default \"%default\"]"),
		make_option(c("--colors_file"),
				help = "Colors file"),
		make_option("--title", default="(unknown title)",
				help="Title for the output figure [default '%default']"),
		make_option("--legend_pos", default="none",
				help="Legend pos [default '%default']")

)

parser <- OptionParser(usage = "stackbar.R [options] environment_file", option_list=option_list,
		description="An interface to create stackbar plots from environment files")

arguments <- parse_args(parser, positional_arguments = TRUE)
options <- arguments$options

if(invalid(options$colors_file)){
	colors_file <- NA
} else {
    colors_file <- options$colors_file
}

# check if the positional argument is set
if(length(arguments$args) != 1) {
	cat("Incorrect number of required positional arguments\n\n")
	print_help(parser)
	stop()
} else {
	environment_file_path <- arguments$args
}

# check if the input file is accessible
if(file.access(environment_file_path) == -1){
	stop(sprintf("Specified file '%s' does not exist", input_file))
}



df <- as.data.frame(read.csv(environment_file_path, header=FALSE, sep="\t"))
names(df) <- c('oligo', 'sample', 'abundance')

df$oligo <- reorder(df$oligo, df$abundance, FUN=sum)
df$oligo <- factor(df$oligo, levels=levels(df$oligo))
df <- df[order(df$oligo), ]

if(!is.na(colors_file)){
    colors_list <- as.data.frame(read.csv(colors_file, header=FALSE, sep="\t"))
    names(colors_list) <- c('oligo', 'color')
    colors_list$oligo <- factor(colors_list$oligo, levels=rev(levels(df$oligo)))
    colors_list <- colors_list[order(colors_list$oligo), ]
} else {
    colors <- rainbow(length(levels(df$oligo)))
    colors <- sample(colors)
    colors_list <- data.frame(cbind(levels(df$oligo), colors))
    names(colors_list) <- c('oligo', 'color')
    print(head(colors_list))
}

P <- function(df){
    p <- ggplot(df, aes(x=factor(sample), y=abundance, fill = oligo))
    p <- p + geom_bar(position="fill", stat = "identity", width=0.90, colour = 'black')
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10))
    p <- p + theme(axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(legend.text = element_text(size = 8, family='mono'), legend.position = options$legend_pos)
    p <- p + labs(x='', y='Oligotypes', title=paste('Distribution of Oligotypes Among Samples for', options$title, sep=' '))
    p <- p + scale_y_continuous(breaks = NULL)
    p <- p + guides(fill = guide_legend(nrow = 25))
    p <- p + theme(plot.title = element_text(hjust=0, vjust=1))
    p <- p + coord_cartesian(ylim=c(-0.01, 1.01))
    p <- p + scale_fill_manual(values = as.vector(colors_list[[2]]))

    print(p)
}

num_samples <- length(levels(factor(df$sample)))
pdf_w <- num_samples / 4
if(num_samples < 32)
    pdf_w <- 8

png_w <- pdf_w * 100
if (png_w > 20000)
    png_w <- 20000

# PDF
pdf_output <- paste(options$output_file_prefix,".pdf",sep="")
pdf(pdf_output, width = pdf_w, height = 8, pointsize = 8, family='Helvetica')
P(df)
sprintf("Stackbar PDF: '%s'", pdf_output)

# PNG
png_output <- paste(options$output_file_prefix,".png",sep="")
png(png_output, width = png_w, height = 800, units = "px", pointsize = 12, bg = "transparent", type = c("cairo", "cairo-png", "Xlib", "quartz"))
P(df)
sprintf("Stackbar PNG: '%s'", png_output)
