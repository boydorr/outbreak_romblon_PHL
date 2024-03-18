# Estimated detection ----

library(data.table) 
library(ggplot2) 
#library(lubridate)
#library(magrittr)
#library(scales)
#library(ggplot2)
#library(patchwork)
library(cowplot)
#library(ggridges)
library(here)
library(ggbeeswarm)
#select <- dplyr::select

# data ----
est_detect <- fread(here("output/tm_trees/detection_ests.csv"))
scenarios <- fread(here("output/tm_trees/scenarios.csv"))

# join up with scenarios ----
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

est_detect <- est_detect[scenarios, on = "scenario"]

# estimated detects ----
est_summ <-
  est_detect[, .(mean_detect = mean(detection)), by = c("scenario", "era", "sim", "data_used", "prune_type", "location")]
est_summ$scenario_label <- paste(est_summ$data_used, est_summ$prune_type,
                                  est_summ$location, sep="\n")
est_summ$scenario_label <- factor(est_summ$scenario_label, 
                                  levels=levels(factor(est_summ$scenario_label))[c(31,32,23,24,27:30,17:22,25,26,
                                                                                 15,16,7,8,11:14,1:6,9,10)])
time_labs <- c("Sep 2022 - Feb 2023", "Mar - Sep 2023", "Overall")
names(time_labs) <- c("First", "Second", "Overall")

est_qs <- # But the 95%CIs are by sim - NOT across sims!
  est_detect[, .(mean_detect = mean(detection),
                 LCI = quantile(detection, 0.025),
                 UCI = quantile(detection, 0.975)), 
             by = c("era", "sim", "data_used", "prune_type", "location")]

ests_all <- 
  ggplot(data = est_summ) +
  geom_density(aes(x = mean_detect, 
                   fill = factor(era), group = factor(era)), 
               alpha = 0.75) +
  scale_fill_brewer(palette = "Dark2", name = "Time period",
                    labels = time_labs) +
  theme_minimal_hgrid(font_size = 12) +
  labs(x = "Estimated detection probability", y = "Density") +
  facet_grid(rows = vars(location, data_used), cols = vars(prune_type),
             labeller = label_wrap_gen(25))
  
ggsave(here("output/figures/detection_all_scenarios.jpeg"), ests_all, 
       height = 12, width = 24, bg = "white")

ests_subset <- 
  ests_all %+% 
  dplyr::filter(est_summ, scenario_label %in% 
                  c("Epi & genetic data\nUnpruned\nBarangay centroids",
                    "Epi & genetic data\nUnpruned\nSimulated locations",
                    "Epi & genetic data\nPruned by time & distance 0.99\nSimulated locations")) %+%
  facet_grid(rows = NULL, cols = vars(scenario_label)) 

ggsave(here("output/figures/detection_subset.jpeg"), ests_subset, 
       height = 4, width = 9, bg = "white")
# 
#   
# # probability of detecting at least one case given a chain size of x
# chain_dets <- expand_grid(est_summ, 
#                           size = seq_len(20))
# chain_dets <- mutate(chain_dets, prob =1 - (1 - mean_detect)^size)
# 
# chain_dets_summary <- 
#   chain_dets %>%
#   group_by(size, era) %>%
#   summarize(mean_prob = mean(prob))
# 
# chains <-
#   ggplot() +
#   geom_segment(data = chain_dets_summary, 
#                aes(x = size - 1, xend = size + 1, 
#                    y = mean_prob, yend = mean_prob, color = factor(era)), 
#                size = 1.25) + 
#   geom_quasirandom(data = chain_dets,
#                    aes(x = size, y = prob, color = factor(era), 
#                        group = era), alpha = 0.25) +
#   labs(aes(x = "Probability of detecting chain of size x")) +
#   scale_color_brewer(palette = "Dark2", name = "Time period") +
#   theme_minimal_hgrid(font_size = 12) +
#   labs(y = "Probability of detecting \nat least one case", x = "Chain size")
# ggsave(here("output/figures/epi_only_barangay_centroids/prob_of_detection_by_chain_size.jpeg"), chains, height = 8, width = 10)

# sfig_detection <- 
#   (recov  / ests_main / chains) + plot_layout(guides = "collect") + 
#   plot_annotation(tag_levels = "A") 
# ggsave("figures/supplement/sfig_detection_main.jpeg", sfig_detection, height = 10, width = 8)

