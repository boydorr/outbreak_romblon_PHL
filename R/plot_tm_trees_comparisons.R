# The code in this script was adapted from the 
# [original in the boydorr/PembaRabies GitHub repository](https://github.com/boydorr/PembaRabies/blob/main/6c.trees_se.R)
# by Carlijn Bogaardt in January 2024.
# The original script is licensed under the MIT license, the full text of which 
# is included below:
#  
#   Copyright (c) 2022 Katie Hampson
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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

ggsave(here("output/figures/figS4_consensus_tree_check.jpeg"), 
       consensus_plots_ts, height = 12, width = 24)

pdf(here("output/figures/figS4_consensus_tree_check.pdf"), width=24, height=12)     
consensus_plots_ts
dev.off()

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

