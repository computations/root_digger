#!/usr/bin/env Rscript

library("data.table")
library("ggplot2")
timeing_data <- fread("./simulated_exp/test_results_times.csv")

trees <- unique(timeing_data$tree)
taxa <- rep("Taxa", length(trees))

aligns <- unique(timeing_data$alignment)
site <- rep("Sites", length(aligns))

timeing_data[, count:=.N, by=list(tree, alignment)]
colMax <- function(data) sapply(data, max, na.rm= TRUE)
data_len = colMax(timeing_data)['count']
timeing_data[, tree_name := paste(tree, "Taxa"), by=list(tree)]
timeing_data[, alignment_name := paste(alignment, "Sites"), by=list(alignment)]

bins = ceiling(2 * data_len^(1/3))

ggplot(data = timeing_data, aes(x = rd_time)) + geom_histogram(aes(y = ..count../sum(..count..)), bins = bins) + 
    facet_wrap(tree_name ~ alignment_name, scales="free") + labs(x = "Time(s)", y = "Density") + theme_minimal()

ggsave("./simulated_exp/rd_time_hist.png")

ggplot(data = timeing_data, aes(x = iq_time)) + geom_histogram(aes(y = ..count../sum(..count..)), bins = bins) + 
    facet_wrap(tree_name ~ alignment_name, scales="free") + labs(x = "Time(s)", y = "Density") + theme_minimal()

ggsave("./simulated_exp/iq_time_hist.png")
