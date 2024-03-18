# Estimating detection probabilities from trees ----

if(Sys.getenv("SLURM_JOB_ID") != "") {
  ncores <- as.numeric(Sys.getenv("SLURM_NTASKS")) 
} else {
  ncores <- parallel::detectCores() - 1
}

cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

# Packages 
library(treerabid)  
library(data.table)
library(foreach)
library(iterators)
library(doRNG)
library(tidyr)
library(here)
library(dplyr)
library(lubridate)

# Estimate detection ----
ttrees <- fread(here("output/tm_trees/trees_all.gz"))
scenarios <- fread(here("output/tm_trees/scenarios.csv"))

# choose 100 random, mcc, & majority for each ----
tree_ids <- ttrees[, .(maj = sum(majority), mcc = sum(mcc)), 
                   by = c("scenario", "sim")][maj > 0 | mcc > 0]
mcc_maj_trees <- ttrees[tree_ids[, .(scenario, sim)], on = c("sim", "scenario")]

# sample trees
set.seed(1677)
selected <- ttrees[!tree_ids[, .(scenario, sim)], on = c("sim", "scenario")][, .(maj = sum(majority)), by = c("sim", "scenario")][, .(sim = sample(sim, 100)), by = "scenario"]
#In the line above I added in a step to reduce to 1 row per sim to avoid duplicate tree-scenario combinations
comp_times_dist <- rbind(mcc_maj_trees, ttrees[selected, on = c("sim", "scenario")])
fwrite(comp_times_dist, here("output/tm_trees/comp_diffs.gz"))

# known and select cols
known <- unique(ttrees[type == "traced"]$id_case)
comp_times_dist[, known_kappa := ifelse(id_case %in% known, 1, 0)]
comp_times_dist <- comp_times_dist[, c("t", "t_diff", "scenario", "sim", "known_kappa")]

# Also want to run analyses for (1) the period with genetic data, and (2) the period without genetic data 
animal_cases <- read.csv(here("data/clean/dat_animal_for_Romblon_tm_trees.csv")) 
last_date <- max(animal_cases[!is.na(animal_cases$samples_sequenced_seq_id),]$date_symptoms_started_proxy)

# Split up and run estimates (overall & by time period)
comp_times_dist[, levs := interaction(sim, 
                                      scenario,
                                      drop = TRUE)]
comp_times_dist <- comp_times_dist[!is.na(t_diff)]
t_diffs <- split(comp_times_dist, comp_times_dist$levs)

est_detect <-
  foreach(k = iter(t_diffs), .combine = rbind) %do% {
    
    split_by_time <- list(overall = k, 
                          first = k[t < ymd(last_date)], 
                          second = k[t > ymd(last_date)])
    ests <- lapply(split_by_time, 
                   function(x) {
                     out <- fit_sims_pi(t_diff = x$t_diff, 
                                        nsims = 10, 
                                        candidate_pis = seq(0.01, 0.99, by = 0.01),
                                        si_fun = treerabid::si_fun_lnorm, 
                                        params = treerabid::params_treerabid, 
                                        alpha = 0.01, 
                                        known_kappas = x$known_kappa,
                                        seed = as.numeric(x$levs[1]))
                     
                     data.table(detection = out, 
                                x[1, -c("known_kappa", "t_diff", "levs")], 
                                nobs = nrow(x))
                   }
    )
    ests[[1]]$era <- "Overall"
    ests[[2]]$era <- "First"
    ests[[3]]$era <- "Second"
    
    rbindlist(ests)
    
  }

parallel::stopCluster(cl)

# Write out files of the trees and the links (consensus & all)
fwrite(est_detect, here("output/tm_trees/detection_ests.csv"))

