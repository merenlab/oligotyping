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
        make_option(c("--otu_limits"), default=F,
                help = "Color 3% OTU limits [default \"%default\"]"),
        make_option(c("-o", "--output_file_prefix"), default="unknown",
                help = "Output file prefix for visualization files [default \"%default\"]"),
        make_option(c("-d", "--distance_col"), default="horn",
                help = "Distance metric for columns [default \"%default\"]"),
        make_option(c("-r", "--distance_row"), default="horn",
                help = "Distance metric for rows [default \"%default\"]"),
        make_option(c("-c", "--clustering"), default="ward",
                help = "Clistering method [default \"%default\"]"),
        make_option(c("--pdf_size"), default=20,
                help = "PDF output height and width [default \"%default\"]"),
        make_option(c("--treeheight_col"), default=100,
                help = "Dendrogram size for columns [default \"%default\"]"),
        make_option(c("--treeheight_row"), default=100,
                help = "Dendrogram size for rows (0 would make it disappear) [default \"%default\"]"),
        make_option("--title", default="(unknown title)",
                help="Title for the output figure [default '%default']")
)

parser <- OptionParser(usage = "heatmap.R [options] input_file", option_list=option_list,
        description="Visualize a distance matrix")

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

drows<-vegdist(raw_data, method='euclidean')
dcols<-vegdist(t(raw_data), method='euclidean', na.rm=TRUE)


P <- function(){
    if(options$otu_limits == TRUE){
        brks1 <- c(0, seq(97, 98.5, length.out = 100))
        cuts1 <- cut(100, breaks = brks1)
        colors1 <- colorRampPalette(c("navy", "#FFFFFF"))(length(levels(cuts1)))

        brks2 <- c(seq(98.500000000001, 100, length.out = 100))
        cuts2 <- cut(100, breaks = brks2)
        colors2 <- colorRampPalette(c("#FFFFFF", "#FF0000"))(length(levels(cuts2)))

        colors <- c("gray", colors1, colors2)
        brks <- c(brks1, brks2)

        pheatmap(raw_data,
                scale="none",
                border_color = NA,
                clustering_distance_rows=drows,
                clustering_distance_cols=dcols,
                clustering_method=options$clustering,
                fontsize_row=2,
                fontsize_col=2,
                color=colors,
                breaks=brks,
                main=options$title, 
                annotation = annotation, 
                treeheight_col = options$treeheight_col, 
                treeheight_row = options$treeheight_row,
                )
    } else {
        pheatmap(raw_data,
                scale="none",
                border_color = NA,
                clustering_distance_rows=drows,
                clustering_distance_cols=dcols,
                clustering_method=options$clustering,
                fontsize_row=2,
                fontsize_col=2,
                main=options$title, 
                annotation = annotation, 
                treeheight_col = options$treeheight_col, 
                treeheight_row = options$treeheight_row,
        )
    }
}

pdf_output <- paste(options$output_file_prefix,".pdf",sep="")
pdf(pdf_output, width = options$pdf_size, height = options$pdf_size, pointsize = 6, family='mono')
P()
sprintf("PDF: '%s'", pdf_output)
dev.off()
