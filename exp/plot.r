#!/usr/bin/env Rscript

library("data.table")
library("ggplot2")
library("cowplot")
library("optparse")

parser <- OptionParser()
parser <- add_option(parser, c("-d", "--data"), type = "character", default = NULL,
    help = "Datafile to make plots with", metavar = "character")

parser <- add_option(parser, c("-p", "--path"), type = "character", default = NULL,
    help = "Directory to place the plots", metavar = "character")

opts = parse_args(parser)

timing_data <- fread(opts$data)

taxa <- unique(timing_data$taxa)

sites <- unique(timing_data$alignment)

timing_data[, `:=`(count, .N), by = list(taxa, sites)]
colMax <- function(data) sapply(data, max, na.rm = TRUE)
timing_data[, `:=`(tree_name, paste(taxa, "Taxa")), by = list(taxa)]
timing_data[, `:=`(alignment_name, paste(sites, "Sites")), by = list(sites)]
timing_data$ordered_tree_name = factor(timing_data$tree_name, levels = c("10 Taxa",
    "20 Taxa", "30 Taxa", "50 Taxa", "100 Taxa"))

sites_selector = "100000 Sites"
tree_selector = "20 Taxa"

summary(timing_data[(taxa=="10" | taxa == "50") & (sites == "1000" | sites == "8000")])
timing_data <- timing_data[(taxa=="10" | taxa == "50") & (sites == "1000" | sites == "8000")]

rd_es_times = timing_data[program == "rd_es"][tree_name == tree_selector][alignment_name ==
    sites_selector]
rd_es_time_median = mean(rd_es_times$time)
rd_no_es_times = timing_data[program == "rd_nes"][tree_name == tree_selector][alignment_name ==
    sites_selector]
rd_no_es_time_median = mean(rd_no_es_times$time)
print(rd_no_es_time_median/rd_es_time_median)

time_boxplot = ggplot(data = timing_data, aes(y = time, x = program, fill = program)) +
    geom_boxplot() + facet_wrap(ordered_tree_name ~ alignment_name, scales = "free",
    shrink = TRUE) + labs(title = "", x = "", y = "Time(s)", fill = "Program") +
    scale_fill_discrete(breaks = c("rd_es", "rd_nes", "iq"), labels = c("RootDigger ES",
        "RootDigger No ES", "IQ-TREE")) + theme_minimal() + theme(axis.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5))

ggsave(paste(opts$path, "/melted_time_boxplot.png", sep=""))

ggplot(data = timing_data, aes(y = root_distance, x = program, fill = program)) +
    geom_boxplot() + facet_grid(ordered_tree_name ~ alignment_name, scales = "free_y") +
    labs(title = "Topological distance of inferred root to true root", x = "", y = "Distance",
        fill = "Program") + scale_fill_discrete(breaks = c("rd_dist", "iq_dist"),
    labels = c("RootDigger", "IQ-TREE")) + theme_minimal() + theme(axis.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5))
ggsave(paste(opts$path, "/melted_dist_box.png", sep=""))

ggplot(data = timing_data, aes(y = normalized_path_distance, x = program, fill = program)) +
    geom_boxplot() + facet_grid(ordered_tree_name ~ alignment_name, scales = "free_y") +
    labs(title = "Path distance of inferred root to true root", x = "", y = "Distance",
        fill = "Program") + scale_fill_discrete(breaks = c("rd_es", "rd_nes", "iq"),
    labels = c("RootDigger ES", "RootDigger No ES", "IQ-TREE")) + theme_minimal() +
    theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5))
ggsave(paste(opts$path, "/melted_path_dist_box.png", sep=""))

distance_boxplot = ggplot(data = timing_data, aes(y = normalized_root_distance, x = program,
    fill = program)) + geom_boxplot() + facet_grid(ordered_tree_name ~ alignment_name,
    scales = "free_y") + labs(title = "", x = "", y = "Distance") + theme_minimal() +
    theme(axis.text.x = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5))
ggsave(paste(opts$path, "/melted_norm_dist_box.png", sep=""))

multi_plot <- plot_grid(distance_boxplot, time_boxplot, labels = "AUTO", rel_widths = c(1,
    1.4))
ggsave(paste(opts$path, "/time_distance_boxplot.png", sep="") , plot = multi_plot, width = 10)
