#!/usr/bin/env Rscript
#
# generates NMDS plots with metadata. here how you run it:
#
# ./script DATA METADATA DISTANCE MAPPING_VARIABLE TITLE OUTPUT_PREFIX
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

library(vegan)
library(gtools)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
csv_path <- args[1]
metadata_path <- args[2]
distance <- args[3]
mapping_variable <- args[4]
title_text <- args[5]
output_file_prefix <- args[6]

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

metadata <- as.data.frame(read.csv(metadata_path, header=TRUE, sep="\t"))
row.names <- metadata$samples
col.names <- colnames(metadata)

mds <- metaMDS(csv[,-1], distance=distance)

NMDS = data.frame(MDS1 = mds$points[,1], MDS2 = mds$points[,2], group=with(metadata, get(mapping_variable)))

NMDS.mean=aggregate(NMDS[,1:2], list(group=with(metadata, get(mapping_variable))), mean)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100){
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
}


df_ell <- data.frame()
only_two_groups <- list()

for(g in levels(NMDS$group)){
    if (nrow(NMDS[NMDS$group==g,]) == 2)
        only_two_groups <- c(only_two_groups, g)
    else
        df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,], veganCovEllipse(cov.wt(cbind(MDS1,MDS2), wt=rep(1/length(MDS1), length(MDS2)))$cov, center=c(mean(MDS1),mean(MDS2))))), group=g))
}

P <- function(){
    if (length(only_two_groups) > 0)
        ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = group)) + geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=0.5, linetype=1) + annotate("text",x=NMDS.mean$MDS1,y=NMDS.mean$MDS2,label=NMDS.mean$group) + geom_line(data = NMDS[NMDS$group %in% only_two_groups,], aes(color = group))
    else
        ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = group)) + geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=0.5, linetype=1) + annotate("text",x=NMDS.mean$MDS1,y=NMDS.mean$MDS2,label=NMDS.mean$group)
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
