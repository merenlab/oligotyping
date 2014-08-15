#!/usr/bin/env Rscript
#
# generates heatmaps.
#

suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(optparse))

# command line options
option_list <- list(
		make_option(c("-m", "--metadata"),
				help = "Metadata file"),
		make_option(c("-o", "--output_file_prefix"), default="unknown",
				help = "Output file prefix for visualization files [default \"%default\"]"),
		make_option(c("-d", "--distance_col"), default="horn",
				help = "Distance metric for columns [default \"%default\"]"),
		make_option(c("-r", "--distance_row"), default="horn",
				help = "Distance metric for rows [default \"%default\"]"),
		make_option(c("-c", "--clustering"), default="ward",
				help = "Clistering method [default \"%default\"]"),
		make_option(c("--pdf_height"), default=9,
				help = "PDF output height [default \"%default\"]"),
		make_option(c("--treeheight_col"), default=100,
				help = "Dendrogram size for columns [default \"%default\"]"),
		make_option(c("--treeheight_row"), default=100,
				help = "Dendrogram size for rows (0 would make it disappear) [default \"%default\"]"),
		make_option(c("--show_rownames"), default=F,
				help = "Show row names [default \"%default\"]"),
		make_option(c("-s", "--scale_the_other_way"), default=F,
				help = "Scale based on columns [default \"%default\"]"),
		make_option("--title", default="(unknown title)",
				help="Title for the output figure [default '%default']")
)

parser <- OptionParser(usage = "heatmap.R [options] input_file", option_list=option_list,
		description="An interface to pheatmap")

arguments <- parse_args(parser, positional_arguments = TRUE)
options <- arguments$options

if(invalid(options$metadata))
	options$metadata <- NA

# check if the positional argument is set
if(length(arguments$args) != 1) {
	cat("Incorrect number of required positional arguments\n\n")
	print_help(parser)
	stop()
} else {
	input_file <- arguments$args
}

# check if the input file is accessible
if(file.access(input_file) == -1){
	stop(sprintf("Specified file '%s' does not exist", input_file))
}

# set annotation
if(is.na(options$metadata)){
	annotation <- NA
} else {
	print(file.access(options$metadata))
	if(file.access(options$metadata) == -1){
		stop(sprintf("Mapping file '%s' does not exist", options$metadata))
	}
	annotation = as.data.frame(read.csv(options$metadata, header=TRUE, row.names = 1, sep="\t"))
}

raw_data <- t(data.matrix(read.table(input_file, header = TRUE, row.names = 1,sep="\t")))

if(options$scale_the_other_way){
    scaled_data <- as.matrix(scale(raw_data), scale = T, center=F)
} else {
    scaled_data <- t(as.matrix(scale(t(raw_data), scale = T, center=F)))
}

drows<-vegdist(raw_data, method=options$distance_row)
dcols<-vegdist(t(raw_data), method=options$distance_col, na.rm=TRUE)

pdf_width <- ncol(scaled_data) / 4
if(pdf_width < 10)
	pdf_width <- 10

if(options$show_rownames == F && nrow(scaled_data) < 70)
	options$show_rownames <- T

P <- function(){
	pheatmap(scaled_data,
			scale="none",
			border_color = NA,
			clustering_distance_rows=drows,
			clustering_distance_cols=dcols,
			clustering_method=options$clustering,
			fontsize_row=8,
			main=options$title, 
			annotation = annotation, 
			treeheight_col = options$treeheight_col, 
			treeheight_row = options$treeheight_row,
			show_rownames = options$show_rownames)
}

pdf_output <- paste(options$output_file_prefix,".pdf",sep="")
pdf(pdf_output, width = pdf_width, height = options$pdf_height, pointsize = 6, family='mono')
P()
sprintf("PDF: '%s'", pdf_output)
dev.off()

png_output <- paste(options$output_file_prefix,".png",sep="")
png(png_output, width = pdf_width * 100, height = options$pdf_height * 100, units = "px", pointsize = 12, bg = "transparent", type = c("cairo", "cairo-png", "Xlib", "quartz"))
P()
sprintf("PNG: '%s'", png_output)
dev.off()


