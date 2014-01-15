#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(gridBase))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(optparse))

# command line options
option_list <- list(
        make_option(c("-o", "--output_file_prefix"), default="unknown",
                help = "Output file prefix [default \"%default\"]"),
        make_option(c("-d", "--distance"), default="canberra",
                help = "Distance metric [default \"%default\"]"),
        make_option(c("-c", "--clustering"), default="complete",
                help = "Clustering method [default \"%default\"]"),
        make_option(c("-O", "--otu_limits"), default=TRUE,
                help = "Mark de facto 97 percent similarity level [default \"%default\"]"),
        make_option(c("-C", "--colors_file"),
                help = "TAB delimited color code for each unit"),
        make_option(c("-v", "--visualize"), default=TRUE,
                help = "Visualize the distance between sequences [default \"%default\"]")
)

parser <- OptionParser(usage = "script.R [options] fasta_file", option_list=option_list,
        description="A script that visualizes similarity between sequences in a FASTA file")

arguments <- parse_args(parser, positional_arguments = TRUE)
options <- arguments$options

# check if the positional argument is set
if(length(arguments$args) != 1) {
    cat("Incorrect number of required positional arguments\n\n")
    print_help(parser)
    stop()
} else {
    input_file_path <- arguments$args[1]
}

if(file.access(input_file_path) == -1){
    stop(sprintf("Input file '%s' does not exist", input_file_path))
}

if(invalid(options$colors_file)){
    colors_file <- NA
} else {
    if(file.access(options$colors_file) == -1){
        stop(sprintf("Colors file '%s' does not exist", options$colors_file))
    }
    colors_file = options$colors_file
}


DIST_BETWEEN_TWO_SEQS <- function(s1, s2){
    #s1 <- 'AAGGCCTACCAAGGCGACGATCAGTAGCGGGTCTGAGAGGATGATCCGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATTTTGGACAATGGGCGCAAGCCTGATCCAGCCATGCCGCGTGTCTGAAGAAGGCCTTCGGGTTGTAAAGGACTTTTGTCAGGGAAGAAAAGGCTGTTGCTAATATCAGCGGCTGATGACGGTACCTGAAGAATAAGCACCGGCTAACTACGTG'
    #s2 <- 'AAGGCCTACCAAGGCGACGATCAGTAGCGGGTCTGAGAGGATGATCCGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATTTTGGACAATGGGCGCAAGCCTGATCCAGCCATGCCGCGTGTCTGAAGAAGGCCTTCGGGTTGTAAAGGACTTTTGTCAGGGAAGAAAAGGCTGTTGCTAATATCGACAGCTGATGACGGTACCTGAAGAATAAGCACCGGCTAACTACGTG'
    #s1 <- 'ATAACG'
    #s2 <- 'ACCG'
    s1 <- gsub("-", "", s1)
    s2 <- gsub("-", "", s2)
    
    alignment <- pairwiseAlignment(s1, s2, gapOpening = -2, gapExtension = -8, scoreOnly = FALSE)
    s1a <- as.character(pattern(alignment))
    s2a <- as.character(subject(alignment))
    mismatch_map <- strsplit(c(s1a, s2a), split= '')
    num_mismatches <- length(which(mismatch_map[[1]] != mismatch_map[[2]]))
    return(100 - (num_mismatches * 100 / nchar(s1a)))
}


OLIGO_DIST <- function(dfx){
    dfx$OLIGO <- factor(dfx$OLIGO)
    
    num_oligos <- nrow(dfx)
    dist_mat <- matrix(nrow=num_oligos, ncol=num_oligos)
    colnames(dist_mat) <- dfx$OLIGO
    rownames(dist_mat) <- dfx$OLIGO

    for(i in 1:num_oligos) {
        for(j in i:num_oligos){
            if(i == j){
                dist_mat[i,j] <- 100.0
                next
            }
            o1 <- dfx$OLIGO[i]
            o2 <- dfx$OLIGO[j]
            d <- DIST_BETWEEN_TWO_SEQS(dfx[dfx$OLIGO == o1, ]$REP_SEQ, dfx[dfx$OLIGO == o2, ]$REP_SEQ)
            dist_mat[i,j] <- d
            dist_mat[j,i] <- d
        }
    }
 
    updated_df <- melt(dist_mat)
    names(updated_df) <- c('OLIGO1', 'OLIGO2', 'DIST')
    
    # find the best order
    d <- vegdist(dist_mat, method=options$distance) 
    fit <- hclust(d, method=options$clustering) #, "single", "complete", "average", "mcquitty", "median" or "centroid"
    oligos_order <- fit$labels[fit$order]
    
    # set the order of samples based on clustering results:
    updated_df$OLIGO1 <- factor(updated_df$OLIGO1, levels=rev(oligos_order))
    updated_df$OLIGO2 <- factor(updated_df$OLIGO2, levels=oligos_order)
    
    dist_mat <- as.matrix(dcast(updated_df, OLIGO1 ~ OLIGO2, value.var="DIST"))
    dist_mat <- apply(dist_mat, 2, rev)
    write.table(dist_mat, file = paste(options$output_file_prefix,".txt",sep=""), sep = "\t", row.names = FALSE, quote = FALSE)

    return(updated_df)
}


#################################################
#####################  MAIN  ####################
#################################################
df <- data.frame(readDNAStringSet(input_file_path, format="fasta"))

if(is.na(colors_file)){
    if(length(row.names(df)) < 64){
        colors <- data.frame(OLIGO=row.names(df), COLOR=sample(colours(), length(row.names(df))), stringsAsFactors=FALSE)
    } else {
        colors <- data.frame(OLIGO=row.names(df), COLOR=rep("#FFFFFF", length(row.names(df))), stringsAsFactors=FALSE)
    }
} else {
    colors <- read.table(colors_file, header = FALSE, sep="\t", comment.char=";")
    names(colors) <- c("OLIGO", "COLOR")
}


dfx <- data.frame(OLIGO=character(),
        REP_SEQ=character(),
        COLOR=character(),
        stringsAsFactors=FALSE)
names(dfx) <- c('OLIGO', 'REP_SEQ', 'COLOR')

N <- 0
for(oligo in row.names(df)){
    N <- N + 1
    dfx[N, ] <- c(OLIGO = oligo, REP_SEQ = as.character(df[oligo, ]), COLOR = as.character(colors[colors$OLIGO == oligo, ]$COLOR))
}

updated_df <- OLIGO_DIST(dfx)

if(options$visualize){
    dfy <- dfx
    dfx$OLIGO <- factor(dfx$OLIGO, levels=levels(updated_df$OLIGO1))
    dfy$OLIGO <- factor(dfy$OLIGO, levels=rev(levels(updated_df$OLIGO1)))
    
    theme_remove_all <- theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_blank(),
                              axis.title.y = element_blank(),
                              axis.text.y = element_blank(),
                              axis.title.x = element_blank(),
                              axis.text.x = element_blank(),
                              axis.ticks = element_blank())
    
    p1 <- ggplot(updated_df, aes(OLIGO1, OLIGO2, fill = DIST))
    p1 <- p1 + geom_tile(width=0.98, height=0.98)
    if(options$otu_limits){
      p1 <- p1 + scale_fill_gradient2(low = "steelblue", high = "red", mid="white", midpoint=98.5, limits=c(97,100), na.value='#AACCAA')
    } else {
      p1 <- p1 + scale_fill_gradient(low = "steelblue", high = "red", na.value='#AACCAA')
    }
    p1 <- p1 + theme_remove_all
    p1 <- p1 + theme(axis.text.x = element_text(angle = 90, size=8, vjust=0, color="black", family="Helvetica"), legend.position = "bottom", panel.background = element_blank())
    p1 <- p1 + xlab(NULL) + ylab(NULL)
    p1 <- p1 + theme(legend.position=c(0.9,0.8))
    p1 <- p1 + theme(legend.background = element_rect())
    p1 <- p1 + theme(legend.background = element_rect(fill=rgb(1,1,1), size=0.5, linetype=2));
    
    colors_dfx <- as.vector(colors[match(levels(dfx$OLIGO), colors$OLIGO), ]$COLOR)
    p2 <- ggplot(dfx, aes(OLIGO, 1, fill=OLIGO))
    p2 <- p2 + geom_tile(width=0.8)
    p2 <- p2 + theme(legend.position = "none", axis.title.y=element_blank(), axis.text.y=element_blank())
    p2 <- p2 + theme_remove_all
    p2 <- p2 + xlab(NULL) + ylab(NULL)
    p2 <- p2 + scale_fill_manual(values = colors_dfx)
    
    colors_dfy <- as.vector(colors[match(levels(dfy$OLIGO), colors$OLIGO), ]$COLOR)
    p3 <- ggplot(dfy, aes(OLIGO, 1, fill=OLIGO))
    p3 <- p3 + geom_tile(width=0.8)
    p3 <- p3 + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank())
    p3 <- p3 + coord_flip()
    p3 <- p3 + theme_remove_all
    p3 <- p3 + xlab(NULL) + ylab(NULL)
    p3 <- p3 + scale_fill_manual(values = colors_dfy)
    
    gp1<- ggplot_gtable(ggplot_build(p1))
    leg <- which(sapply(gp1$grobs, function(x) x$name) == "guide-box")
    legend <- gp1$grobs[[leg]]
    
    gp1<- ggplot_gtable(ggplot_build(p1))
    gp2<- ggplot_gtable(ggplot_build(p2))
    gp3<- ggplot_gtable(ggplot_build(p3))
    
    maxWidth = unit.pmax(gp1$widths[2:3], gp2$widths[2:3])
    maxHeight <- unit.pmax(gp1$heights[4:5], gp3$heights[4:5])
    gp1$widths[2:3] <- maxWidth
    gp2$widths[2:3] <- maxWidth
    gp1$heights[4:5] <- maxHeight
    gp3$heights[4:5] <- maxHeight
    
    blank <- grid.rect(gp=gpar(col="white"))
    pdf_output <- paste(options$output_file_prefix,".pdf",sep="")
    pdf(pdf_output, width=10, height=10)
    grid.arrange(blank, gp2, gp3, gp1, heights=c(1/10, 9/10), widths=c(1/10, 9/10), ncol=2)
    dev.off()
    
    png_output <- paste(options$output_file_prefix,".png",sep="")
    png(png_output, width=1200, height=1200)
    grid.arrange(blank, gp2, gp3, gp1, heights=c(1/10, 9/10), widths=c(1/10, 9/10), ncol=2)
    dev.off()
}
