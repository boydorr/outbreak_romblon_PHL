# Comparison of consensus trees across 28 scenarios ----

library(ggplot2)
library(lubridate)
library(dplyr)
library(ggraph)
library(igraph)
library(patchwork)
library(sf)
library(ggthemes)
library(ggforce)
library(here)
library(data.table) 
source(here("R/plot_lineage_ts.R"))

# Loading data
tree_consensus <- fread(here("output/tm_trees/links_consensus_consistent.csv"))
scenarios <- fread(here("output/tm_trees/scenarios.csv"))
animal_cases <- fread(here("data/clean/dat_animal_for_Romblon_tm_trees.csv"))
case_dates <- data.table(id_case = animal_cases$id,
                        symptoms_started = animal_cases$date_symptoms_started_proxy)

# Cleaning up scenario labels
scenarios[, c("data_used", "prune_type", "location") := .(ifelse(!use_gen, "Epi data only", "Epi & genetic data"), 
                                                          fcase(!prune, "Unpruned", 
                                                                time_cutoff == 0.99  & dist_cutoff == 0.99, 
                                                                "Pruned by time & distance 0.99", 
                                                                time_cutoff == 0.975 & dist_cutoff == 0.975, 
                                                                "Pruned by time & distance 0.975", 
                                                                time_cutoff == 0.95 & dist_cutoff == 0.95, 
                                                                "Pruned by time & distance 0.95", 
                                                                time_cutoff == 0.975 & dist_cutoff == 0.99, 
                                                                "Pruned by time 0.975 & distance 0.99", 
                                                                time_cutoff == 0.99 & dist_cutoff == 1, 
                                                                "Pruned by time 0.99",
                                                                time_cutoff == 0.975 & dist_cutoff == 1, 
                                                                "Pruned by time 0.975",
                                                                time_cutoff == 0.95 & dist_cutoff == 1, 
                                                                "Pruned by time 0.95"),
                                                          ifelse(location == "centroid", "Barangay centroids", "Simulated locations"))]
scenarios <- scenarios[, c("scenario", "data_used", "prune_type", "use_gen", "location")]
scenarios$prune_type <- factor(scenarios$prune_type, 
                               levels = c("Unpruned", "Pruned by time 0.99",
                                          "Pruned by time 0.975", "Pruned by time 0.95",
                                          "Pruned by time & distance 0.99",
                                          "Pruned by time & distance 0.975", 
                                          "Pruned by time & distance 0.95",
                                          "Pruned by time 0.975 & distance 0.99"))
scenarios$data_used <- factor(scenarios$data_used, 
                              levels = c("Epi data only", "Epi & genetic data"))
scenarios$location <- factor(scenarios$location)

# Plotting
tree_consensus <- tree_consensus[scenarios, on = "scenario"]
con_plots <- 
  lapply(split(tree_consensus,
               interaction(tree_consensus$prune_type, 
                           tree_consensus$data_used,
                           tree_consensus$location)), 
         function(x) {
           if(x$data_used[1] == "Epi data only") {
             x <- x[, -c("lineage_chain", "lineage")]
             x <- x[animal_cases[, c("id","lineage")], on = c(id_case="id")]
             x$lineage_chain <- x$lineage
           } else {colour_var = "lineage_chain"}
           out <-  plot_lin_ts(x, case_dates,
                               colour_by = "lineage_chain",
                               colour_edges = FALSE,
                               shape_unsampled = 21,
                               id_labels = NA,
                               timeline_breaks = 3)
           out <- out + labs(subtitle = paste0(x$data_used[1], "\n", 
                                               x$prune_type, "\n", x$location))
         }
  )

consensus_plots_ts <- 
  wrap_plots(con_plots, ncol = 8, nrow = 4) + 
  plot_annotation(title = "Consensus tree comparison") &
  panel_border(size = 0.1, color = "black")

ggsave(here("output/figures/consensus_tree_check.jpeg"), 
       consensus_plots_ts, height = 12, width = 24)

# Comparison of three selected scenario trees ----

selected_plots <-
  con_plots[c("Unpruned.Epi & genetic data.Barangay centroids",
            "Unpruned.Epi & genetic data.Simulated locations",
            "Pruned by time & distance 0.99.Epi & genetic data.Simulated locations")]

consensus_plots_subset <- 
  wrap_plots(selected_plots, ncol = 3, nrow = 1) + 
  plot_annotation(title = "Consensus tree comparison") &
  panel_border(size = 0.1, color = "black")

ggsave(here("output/figures/consensus_tree_subset.jpeg"), 
       consensus_plots_subset, height = 4, width = 9)

