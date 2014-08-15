#!/usr/bin/env Rscript
#
# generates NMDS plots with metadata.
#
# DATA is a TAB delimited contingency table, which should look like this:
#
#   samples     species1     species2       species3        ...
#   sample1        0.4          1.6            7.7          ...
#   sample2        9.9          1.1            9.5          ...
#   sample3        9.1          6.7            5.5          ...
#   ...
#
#   METADATA is the TAB delimited mapping file, like this one:
#
#   samples     category1     category2     ...
#   sample1       deep          saline      ...
#   sample2       deep          fresh       ...
#   sample3      shallow        fresh       ...  
#
#
#   DISTANCE is the distance function (one of these: "manhattan", "euclidean",
#   "canberra", "bray", "kulczynski", "jaccard", "gower", "morisita", "horn", 
#   "mountford", "raup" , "binomial" or "chao")
#
#   MAPPING_VARIABLE one of the categories from the METADATA file.
#

suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))

# command line options
option_list <- list(
		make_option(c("-o", "--output_file_prefix"), default="unknown",
				help = "Output file prefix [default \"%default\"]"),
		make_option(c("-d", "--distance"), default="horn",
				help = "Distance metric [default \"%default\"]"),
		make_option(c("-m", "--mapping_variable"),
				help = "Column in the metadata for sample mapping"),
		make_option("--title", default="(unknown title)",
				help="Title for the output figure [default '%default']")
)

parser <- OptionParser(usage = "script.R [options] input_matrix metadata", option_list=option_list,
		description="A script to generate MDS plots with sample mapping")

arguments <- parse_args(parser, positional_arguments = TRUE)
options <- arguments$options

# check if the positional argument is set
if(length(arguments$args) != 2) {
	cat("Incorrect number of required positional arguments\n\n")
	print_help(parser)
	stop()
} else {
	input_file_path <- arguments$args[1]
	metadata_path <- arguments$args[2]
}

if(file.access(input_file_path) == -1){
	stop(sprintf("Input file '%s' does not exist", input_file_path))
}

if(file.access(metadata_path) == -1){
	stop(sprintf("Metadata file '%s' does not exist", metadata_path))
}

if(invalid(options$mapping_variable))
	stop(sprintf("You must define a mapping variable (-m)"))

data <- as.data.frame(read.table(input_file_path, header = TRUE, sep="\t", comment.char = '&'))
metadata <- as.data.frame(read.table(metadata_path, header=TRUE, sep="\t", comment.char = '&'))

if(names(data)[1] != 'samples')
	stop(sprintf("Data file '%s' does not seem to be formatted properly", input_file_path))
if(names(metadata)[1] != 'samples')
	stop(sprintf("Metadata file '%s' does not seem to be formatted properly", metadata_path))

samples_in_both <- intersect(data$samples, metadata$samples)
data <- data[data$samples %in% samples_in_both, ]
metadata <- metadata[metadata$samples %in% samples_in_both, ]

if(dim(metadata)[1] == 0)
	stop("Samples in the input matrix and metadata do not seem to correspond")

metadata <- metadata[match(data$samples, metadata$samples),]
data$samples <- factor(data$samples)
metadata$samples <- factor(metadata$samples)

if(!options$mapping_variable %in% names(metadata)){
	stop(sprintf("Metadata file does not contain mapping variable '%s'", options$mapping_variable))
}

mds <- metaMDS(data[,-1], distance=options$distance)

NMDS = data.frame(MDS1 = mds$points[,1], MDS2 = mds$points[,2], group=with(metadata, get(options$mapping_variable)))

NMDS.mean=aggregate(NMDS[,1:2], list(group=with(metadata, get(options$mapping_variable))), mean)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100){
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
}


df_ell <- data.frame()
only_two_groups <- list()

for(g in levels(NMDS$group)){
    if (nrow(NMDS[NMDS$group==g,]) < 3)
        only_two_groups <- c(only_two_groups, g)
    else
        df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,], veganCovEllipse(cov.wt(cbind(MDS1,MDS2), wt=rep(1/length(MDS1), length(MDS2)))$cov, center=c(mean(MDS1),mean(MDS2))))), group=g))
}

P <- function(){
    if (length(only_two_groups) > 0){
        p <- ggplot(data = NMDS, aes(MDS1, MDS2))
		p <- p + geom_point(aes(color = group))
		p <- p + geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=0.5, linetype=1)
		p <- p + geom_line(data = NMDS[NMDS$group %in% only_two_groups,], aes(color = group))
		p <- p + annotate("text",x=NMDS.mean$MDS1,y=NMDS.mean$MDS2,label=NMDS.mean$group, size=8)
		p <- p + ggtitle(options$title)
		p <- p + theme_bw()
	}
    else{
		p <- ggplot(data = NMDS, aes(MDS1, MDS2))
		p <- p + geom_point(aes(color = group))
		p <- p + geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=0.5, linetype=1)
		p <- p + annotate("text",x=NMDS.mean$MDS1,y=NMDS.mean$MDS2,label=NMDS.mean$group, size=8)
		p <- p + ggtitle(options$title)
		p <- p + theme_bw()
	}
	
	print(p)
}

# PDF
pdf_output <- paste(options$output_file_prefix,".pdf",sep="")
pdf(pdf_output, width = 16, height = 10, family='Helvetica')
P()
sprintf("Clustering result PDF: '%s'", pdf_output)
dev.off()

# PNG
png_output <- paste(options$output_file_prefix,".png",sep="")
png(png_output, width = 800, height = 600, units = "px", pointsize = 12, bg = "transparent", type = c("cairo", "cairo-png", "Xlib", "quartz"))
P()
sprintf("Clustering result PNG: '%s'", png_output)
dev.off()
