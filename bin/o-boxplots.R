#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(compute.es))

# command line options
option_list <- list(
		make_option(c("--mapping_var"),
				help = "Mapping variable to use from the sample mapping file. If none specified, the first column is used."),
		make_option(c("--output_directory"), default=".",
				help = "Output directory to store images [default \"%default\"]"),
		make_option(c("--remove_outliers"), default=F,
				help = "Remove upper outliers for better scaling [default \"%default\"]"),
		make_option(c("--pdf_height"), default=9,
				help = "PDF output height [default \"%default\"]"))


parser <- OptionParser(usage = "this_script.R [options] input_matrix sample_mapping", option_list=option_list,
		description="Visualize the distribution of a unit among mapping categories")

arguments <- parse_args(parser, positional_arguments = TRUE)
options <- arguments$options

# check if the positional argument is set
if(length(arguments$args) != 2) {
	cat("Incorrect number of required positional arguments\n\n")
	print_help(parser)
	stop()
}

input_matrix <- arguments$args[1]
sample_mapping <- arguments$args[2]

# check if the input file is accessible
if(file.access(input_matrix) == -1){
	stop(sprintf("Input matrix '%s' does not exist", input_matrix))
}
if(file.access(sample_mapping) == -1){
	stop(sprintf("Sample mapping '%s' does not exist", sample_mapping))
}

data <- as.data.frame(read.csv(input_matrix, header=TRUE, sep="\t"))
metadata <- as.data.frame(read.csv(sample_mapping, header=TRUE, sep="\t"))

metadata <- subset(metadata, metadata$samples %in% data$samples)
metadata$samples <- factor(metadata$samples)
metadata <- metadata[with(metadata, order(samples)),] 
data <- data[with(data, order(samples)),] 

if(invalid(options$mapping_var)){
	mapping_variable <- names(metadata)[2]
} else {
	if(options$mapping_var %in% names(metadata)){
		mapping_variable <- options$mapping_var
	} else {
		stop(sprintf("Mapping variable '%s' does not exist in the mapping file", options$mapping_var))
	}
}

pdf_width <- length(levels(metadata[[mapping_variable]]))
if(pdf_width < 4)
	pdf_width <- 4

nodes <- names(data[, !names(data) %in% c("samples")]) 

sums <- rep(NA, length(nodes))
for(i in 1:length(nodes)){
	sums[[i]] <- c(sum(data[[nodes[i]]]))
}
ord <- order(sums, decreasing = TRUE)

get_first_two_factors_and_effect_size <- function(d, mapping_variable, node){
	#Êperform lda on the node with respect to classes
	lda_result <- suppressWarnings(lda(d[[mapping_variable]] ~ ., data = data.frame(d[[node]]), tol=0.00001))
	# find the two factors with the highes mean for this node
	factors_means_df <- data.frame(levels(factor(d[[mapping_variable]])), as.numeric(lda_result$means))
	colnames(factors_means_df) <- c('factors', 'means')
	ordered_factors_means_df <- factors_means_df[with(factors_means_df, order(-means)), ]
	factor_1 <- as.character(ordered_factors_means_df[1,]$factors)
	factor_2 <- as.character(ordered_factors_means_df[2,]$factors)
	
	#Êlda scaling
	w <- lda_result$scaling[,1]
	scaled_frequencies <- as.matrix(d[[node]])%*%w
	effect_size <- abs(mean(scaled_frequencies[d[[mapping_variable]]==factor_1]) - mean(scaled_frequencies[d[[mapping_variable]]==factor_2]))
	return(data.frame(factor_1, factor_2, effect_size))
}

mapping <- with(metadata, get(mapping_variable))
mapping <- factor(mapping, levels=levels(mapping))

for (i in 1:length(nodes[ord])){
	node <- nodes[ord][i]
	
	#Êget the subset of data with the mapping and node of interest i = 1
	d <- data.frame(data[[node]], metadata[[mapping_variable]])
	names(d) <- c(node, mapping_variable)
	mapping <- d[[mapping_variable]]
	mapping <- factor(mapping, levels=levels(mapping))
	
	if(options$remove_outliers == T){
		upper_whisker <- boxplot.stats(d[[node]])$stats[c(1, 5)][2] * 2
		outlier_limit <- mean(d[d[[node]] > upper_whisker, ][[node]])
		if(!is.nan(outlier_limit)){
			d <- d[d[[node]] < outlier_limit, ]
			mapping <- d[[mapping_variable]]
			mapping <- factor(mapping, levels=levels(mapping))
		}
	}
	
	k <- kruskal.test(d[[node]], g = mapping)
	kruskal_p <- k$p.value
	kruskal_s <- as.numeric(k$statistic)
	
	#if (kruskal_p > 0.05)
	#	next
	
	lda_result <- get_first_two_factors_and_effect_size(d, mapping_variable, node)
	factor_1 <- as.character(lda_result$factor_1)
	factor_2 <- as.character(lda_result$factor_2)
	lda_effect_size <- lda_result$effect_size
	
	#Êfactor 1 is the highest mean. lets proceed with that.
	d.test <- d[which(d[[mapping_variable]] %in% factor_1), ][[node]]
	d.ctrl <- d[which(d[[mapping_variable]] %in% factor_2), ][[node]]
	d.mes_results <- mes(mean(d.test), mean(d.ctrl), sd(d.test), sd(d.ctrl), length(d.test), length(d.ctrl))
	
	mes_d <- as.numeric(d.mes_results$MeanDifference[1])
	mes_d_var <- as.numeric(d.mes_results$MeanDifference[2])
	mes_g <- as.numeric(d.mes_results$MeanDifference[3])
	mes_g_var <- as.numeric(d.mes_results$MeanDifference[4])
	
	
	w <- wilcox.test(d.test, d.ctrl)
	wilcox_p <- w$p.value
	wilcox_s <- w$statistic
	
	#if (wilcox_p > 0.05)
	#	next
		
	out_text <- paste(i, node, factor_1, signif(kruskal_s, digits = 3), signif(kruskal_p, digits = 3), signif(wilcox_p, digits = 3), signif(lda_effect_size, digits = 3), sep="\t")
	cat(paste(out_text, '\n', sep=''))
	label <- paste("Kruskal test\np: ", signif(kruskal_p, digits = 2), ', Chi-S: ', signif(kruskal_s, digits = 3), sep = "")	
	label <- paste(label, "\n\nWilcoxon test\np: ", signif(wilcox_p, digits = 2), ', W: ', signif(wilcox_s, digits = 3), sep = "")	
	label <- paste(label, "\n\nMeans To Effect Size -> ", factor_1, "\n", sep="")
	label <- paste(label, "d: ", signif(mes_d, digits = 2), ", var(d): ", signif(mes_d_var, digits = 2), sep="")
	
	g <- ggplot(data=d, aes(mapping, d[[node]], color=mapping))
	g <- g + geom_boxplot() 
	g <- g + geom_jitter(position = position_jitter(width = .05, height = 0), size=4, alpha=0.7) 
	g <- g + theme(axis.text.x = element_text(angle = 90, size=12, colour='black')) 
	g <- g + theme(legend.position = 'none')
	g <- g + labs(x='', y='Percent Abundance', title = node)
	#g <- g + annotate("text", x=(length(levels(factor(mapping))) + 1) / 2, y=max(d[[node]]) * 0.8, label=label, size=4)
	
	pdf_output <- paste(options$output_directory, '/', sprintf("%03d", i), '_', node, ".pdf",sep="")
	pdf(pdf_output, width = pdf_width, height = 8, pointsize = 8)#, family='mono')
	print(g)
	dev.off()
	sprintf("PDF: '%s'", pdf_output)
}
