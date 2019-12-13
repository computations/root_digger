#!/usr/bin/env Rscript

library("data.table")
library("ggplot2")
json_data <- fread("./simulated_exp/test_results_times.csv")

trees <- unique(json_data$tree)
taxa <- rep("Taxa", length(trees))

aligns <- unique(json_data$alignment)
site <- rep("Sites", length(aligns))

json_data[, count:=.N, by=list(tree, alignment)]
colMax <- function(data) sapply(data, max, na.rm= TRUE)
data_len = colMax(json_data)['count']
json_data[, tree_name := paste(tree, "Taxa"), by=list(tree)]
json_data[, alignment_name := paste(alignment, "Sites"), by=list(alignment)]

bins = ceiling(2 * data_len^(1/3))

ggplot(data = json_data, aes(x = rd_time)) + geom_histogram(aes(y = ..count../sum(..count..)), bins = bins) + 
    facet_grid(tree_name ~ alignment_name) + labs(x = "Time(s)", y = "Density")

ggsave("ggtest.png")
