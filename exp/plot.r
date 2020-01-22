#!/usr/bin/env Rscript

library("data.table")
library("ggplot2")
library("cowplot")
timing_data <- fread("./simulated_exp/test_results.csv")

taxa <- unique(timing_data$taxa)

sites <- unique(timing_data$alignment)

timing_data[, `:=`(count, .N), by = list(taxa, sites)]
colMax <- function(data) sapply(data, max, na.rm = TRUE)
timing_data[, `:=`(tree_name, paste(taxa, "Taxa")), by = list(taxa)]
timing_data[, `:=`(alignment_name, paste(sites, "Sites")), by = list(sites)]

time_boxplot = ggplot(data = timing_data, aes(y = time, x = program, fill = program)) + 
    geom_boxplot() + facet_wrap(tree_name ~ alignment_name, scales = "free", shrink = TRUE) + 
    labs(title = "", x = "", y = "Time(s)", fill = "Program") + scale_fill_discrete(breaks = c("rd_es", 
    "rd_nes", "iq"), labels = c("RootDigger ES", "RootDigger No ES", "IQ-TREE")) + 
    theme_minimal() + theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5))

ggsave("./simulated_exp/melted_time_boxplot.png")

ggplot(data = timing_data, aes(y = root_distance, x = program, fill = program)) + 
    geom_boxplot() + facet_grid(tree_name ~ alignment_name, scales = "free_y") + 
    labs(title = "Topological distance of inferred root to true root", x = "", y = "Distance", 
        fill = "Program") + scale_fill_discrete(breaks = c("rd_dist", "iq_dist"), 
    labels = c("RootDigger", "IQ-TREE")) + theme_minimal() + theme(axis.text.x = element_blank(), 
    plot.title = element_text(hjust = 0.5))
ggsave("./simulated_exp/melted_dist_box.png")

ggplot(data = timing_data, aes(y = normalized_path_distance, x = program, fill = program)) + 
    geom_boxplot() + facet_grid(tree_name ~ alignment_name, scales = "free_y") + 
    labs(title = "Path distance of inferred root to true root", x = "", y = "Distance", 
        fill = "Program") + scale_fill_discrete(breaks = c("rd_es", "rd_nes", "iq"), 
    labels = c("RootDigger ES", "RootDigger No ES", "IQ-TREE")) + theme_minimal() + 
    theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5))
ggsave("./simulated_exp/melted_path_dist_box.png")

distance_boxplot = ggplot(data = timing_data, aes(y = normalized_root_distance, x = program, 
    fill = program)) + geom_boxplot() + facet_grid(tree_name ~ alignment_name, 
    scales = "free_y") + labs(title = "", x = "", y = "Distance") + theme_minimal() + 
    theme(axis.text.x = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5))
ggsave("./simulated_exp/melted_norm_dist_box.png")

multi_plot <- plot_grid(distance_boxplot, time_boxplot, labels = "AUTO", rel_widths = c(1, 
    1.4))
ggsave("./simulated_exp/time_distance_boxplot.png", plot = multi_plot, width = 10)
