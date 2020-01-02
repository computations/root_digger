#!/usr/bin/env Rscript

library("data.table")
library("ggplot2")
timing_data <- fread("./simulated_exp/test_results_times.csv")

trees <- unique(timing_data$tree)
taxa <- rep("Taxa", length(trees))

aligns <- unique(timing_data$alignment)
site <- rep("Sites", length(aligns))

timing_data[, count:=.N, by=list(tree, alignment)]
colMax <- function(data) sapply(data, max, na.rm= TRUE)
data_len = colMax(timing_data)['count']
timing_data[, tree_name := paste(tree, "Taxa"), by=list(tree)]
timing_data[, alignment_name := paste(alignment, "Sites"), by=list(alignment)]

bins = ceiling(2 * data_len^(1/3))

ggplot(data = timing_data, aes(x = rd_time)) + geom_histogram(aes(y = ..count../sum(..count..)), bins = bins) + 
    facet_wrap(tree_name ~ alignment_name, scales="free") + labs(x = "Time(s)", y = "Density") + theme_minimal()

ggsave("./simulated_exp/rd_time_hist.png")

ggplot(data = timing_data, aes(x = iq_time)) + geom_histogram(aes(y = ..count../sum(..count..)), bins = bins) + 
    facet_wrap(tree_name ~ alignment_name, scales="free") + labs(x = "Time(s)", y = "Density") + theme_minimal()

ggsave("./simulated_exp/iq_time_hist.png")

distance_data <- fread("./simulated_exp/test_results_dists.csv")

trees <- unique(distance_data$tree)
taxa <- rep("Taxa", length(trees))

aligns <- unique(distance_data$alignment)
site <- rep("Sites", length(aligns))

distance_data[, count:=.N, by=list(tree, alignment)]
data_len = colMax(distance_data)['count']
distance_data[, tree_name := paste(tree, "Taxa"), by=list(tree)]
distance_data[, alignment_name := paste(alignment, "Sites"), by=list(alignment)]
bins = ceiling(2 * data_len^(1/3))

ggplot(data = distance_data, aes(x = rd_dist)) + geom_histogram(aes(y = ..count../sum(..count..)), bins = bins) + 
    facet_wrap(tree_name ~ alignment_name, scales="free") + labs(x = "Dist", y = "Density") + theme_minimal()

ggsave("./simulated_exp/rd_dist_hist.png")

ggplot(data = distance_data, aes(x = iq_dist)) + geom_histogram(aes(y = ..count../sum(..count..)), bins = bins) + 
    facet_wrap(tree_name ~ alignment_name, scales="free") + labs(x = "Dist", y = "Density") + theme_minimal()

ggsave("./simulated_exp/iq_dist_hist.png")

ggplot(data = distance_data, aes(x = rd_path_dist)) + geom_histogram(aes(y = ..count../sum(..count..)), bins = bins) + 
    facet_wrap(tree_name ~ alignment_name, scales="free") + labs(x = "Dist", y = "Density") + theme_minimal()

ggsave("./simulated_exp/rd_path_dist_hist.png")

ggplot(data = distance_data, aes(x = iq_path_dist)) + geom_histogram(aes(y = ..count../sum(..count..)), bins = bins) + 
    facet_wrap(tree_name ~ alignment_name, scales="free") + labs(x = "Dist", y = "Density") + theme_minimal()

ggsave("./simulated_exp/iq_path_dist_hist.png")
