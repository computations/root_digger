#!/usr/bin/env Rscript

library("data.table")
library("ggplot2")
library("cowplot")
timing_data <- fread("./simulated_exp/test_results_times.csv")

trees <- unique(timing_data$tree)
taxa <- rep("Taxa", length(trees))

aligns <- unique(timing_data$alignment)
site <- rep("Sites", length(aligns))

timing_data[, `:=`(count, .N), by = list(tree, alignment)]
colMax <- function(data) sapply(data, max, na.rm = TRUE)
data_len = colMax(timing_data)["count"]
timing_data[, `:=`(tree_name, paste(tree, "Taxa")), by = list(tree)]
timing_data[, `:=`(alignment_name, paste(alignment, "Sites")), by = list(alignment)]

bins = ceiling(2 * data_len^(1/3))

ggplot(data = timing_data, aes(x = rd_time)) + geom_histogram(aes(y = ..count../sum(..count..)), 
    bins = bins) + facet_wrap(tree_name ~ alignment_name, scales = "free") + labs(x = "Time(s)", 
    y = "Density") + theme_minimal()

ggsave("./simulated_exp/rd_time_hist.png")

ggplot(data = timing_data, aes(x = iq_time)) + geom_histogram(aes(y = ..count../sum(..count..)), 
    bins = bins) + facet_wrap(tree_name ~ alignment_name, scales = "free") + labs(x = "Time(s)", 
    y = "Density") + theme_minimal()

ggsave("./simulated_exp/iq_time_hist.png")

melted_time_data = melt(timing_data, measure.vars = c("rd_time", "iq_time"), id.vars = c("tree", 
    "alignment", "tree_name", "alignment_name"))
melted_time_data$ordered_tree_names = factor(melted_time_data$tree_name, levels = c("10 Taxa", 
    "50 Taxa", "100 Taxa"))

time_boxplot = ggplot(data = melted_time_data, aes(y = value, x = variable, fill = variable)) + 
    geom_boxplot() + facet_wrap(ordered_tree_names ~ alignment_name, scales = "free", shrink=TRUE) + 
    labs(title = "", x = "", y = "Time(s)", fill = "Program") + scale_fill_discrete(breaks = c("rd_time", 
    "iq_time"), labels = c("RootDigger", "IQ-TREE")) + theme_minimal() + theme(axis.text.x = element_blank(), 
    plot.title = element_text(hjust = 0.5))

ggsave("./simulated_exp/melted_time_boxplot.png")

distance_data <- fread("./simulated_exp/test_results_dists.csv")
distance_data <- distance_data[order(tree, alignment)]
distance_data$rd_normed_distance = distance_data$rd_dist/distance_data$tree
distance_data$iq_normed_distance = distance_data$iq_dist/distance_data$tree

trees <- unique(distance_data$tree)
taxa <- rep("Taxa", length(trees))

aligns <- unique(distance_data$alignment)
site <- rep("Sites", length(aligns))

distance_data[, `:=`(count, .N), by = list(tree, alignment)]
data_len = colMax(distance_data)["count"]
distance_data[, `:=`(tree_name, paste(tree, "Taxa")), by = list(tree)]
distance_data[, `:=`(alignment_name, paste(alignment, "Sites")), by = list(alignment)]
bins = ceiling(2 * data_len^(1/3))

ggplot(data = distance_data, aes(x = rd_dist)) + geom_histogram(aes(y = ..count../sum(..count..)), 
    bins = bins) + facet_wrap(tree_name ~ alignment_name, scales = "free") + labs(x = "Dist", 
    y = "Density") + theme_minimal()

ggsave("./simulated_exp/rd_dist_hist.png")

ggplot(data = distance_data, aes(x = iq_dist)) + geom_histogram(aes(y = ..count../sum(..count..)), 
    bins = bins) + facet_wrap(tree_name ~ alignment_name, scales = "free") + labs(x = "Dist", 
    y = "Density") + theme_minimal()

ggsave("./simulated_exp/iq_dist_hist.png")

ggplot(data = distance_data, aes(x = rd_path_dist)) + geom_histogram(aes(y = ..count../sum(..count..)), 
    bins = bins) + facet_wrap(tree_name ~ alignment_name, scales = "free") + labs(x = "Dist", 
    y = "Density") + theme_minimal()

ggsave("./simulated_exp/rd_path_dist_hist.png")

ggplot(data = distance_data, aes(x = iq_path_dist)) + geom_histogram(aes(y = ..count../sum(..count..)), 
    bins = bins) + facet_wrap(tree_name ~ alignment_name, scales = "free") + labs(x = "Dist", 
    y = "Density") + theme_minimal()

ggsave("./simulated_exp/iq_path_dist_hist.png")

melted_distance_data = melt(distance_data, measure.vars = c("rd_dist", "iq_dist"), 
    id.vars = c("tree", "alignment", "tree_name", "alignment_name"))
melted_path_distance_data = melt(distance_data, measure.vars = c("rd_path_dist", 
    "iq_path_dist"), id.vars = c("tree", "alignment", "tree_name", "alignment_name"))
melted_normed_distance_data = melt(distance_data, measure.vars = c("rd_normed_distance", 
    "iq_normed_distance"), id.vars = c("tree", "alignment", "tree_name", "alignment_name"))

melted_normed_distance_data$ordered_tree_names = factor(melted_normed_distance_data$tree_name, 
    levels = c("10 Taxa", "50 Taxa", "100 Taxa"))

variable_names <- list(rd_dist = "RootDigger Distance", iq_dist = "IQ-TREE Distance", 
    rd_path_dist = "RootDigger Path Distance", iq_path_dist = "IQ-TREE Path Distance", 
    rd_normed_distance = "RootDigger Normalized Distance", iq_normed_distance = "IQ-TREE Noramlized Distance")

v_l <- function(variable, value) {
    return(variable_names[value])
}

ggplot(data = melted_distance_data, aes(y = value, x = variable, fill = variable)) + 
    geom_boxplot() + facet_grid(tree_name ~ alignment_name, scales = "free_y") + 
    labs(title = "Topological distance of inferred root to true root", x = "", y = "Distance", 
        fill = "Program") + scale_fill_discrete(breaks = c("rd_dist", "iq_dist"), 
    labels = c("RootDigger", "IQ-TREE")) + theme_minimal() + theme(axis.text.x = element_blank(), 
    plot.title = element_text(hjust = 0.5))
ggsave("./simulated_exp/melted_dist_box.png")

ggplot(data = melted_path_distance_data, aes(y = value, x = variable, fill = variable)) + 
    geom_boxplot() + facet_grid(tree_name ~ alignment_name, scales = "free_y") + 
    labs(title = "Path distance of inferred root to true root", x = "", y = "Distance", 
        fill = "Program") + scale_fill_discrete(breaks = c("rd_path_dist", "iq_path_dist"), 
    labels = c("RootDigger", "IQ-TREE")) + theme_minimal() + theme(axis.text.x = element_blank(), 
    plot.title = element_text(hjust = 0.5))
ggsave("./simulated_exp/melted_path_dist_box.png")

distance_boxplot = ggplot(data = melted_normed_distance_data, aes(y = value, x = variable, 
    fill = variable)) + geom_boxplot() + facet_grid(ordered_tree_names ~ alignment_name, scales = "free_y") + 
    labs(title = "", x = "", y = "Distance") + theme_minimal() + theme(axis.text.x = element_blank(), 
    legend.position = "none", plot.title = element_text(hjust = 0.5))
ggsave("./simulated_exp/melted_norm_dist_box.png")

multi_plot <- plot_grid(distance_boxplot, time_boxplot, labels = "AUTO", rel_widths = c(1, 
    1.4))
ggsave("./simulated_exp/time_distance_boxplot.png", plot = multi_plot, width = 10)
