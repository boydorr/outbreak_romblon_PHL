# The code in this script was adapted from the 
# [original in the boydorr/PembaRabies GitHub repository](https://github.com/boydorr/PembaRabies/blob/main/6.TransmissionTrees.R)
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

library(here)
library(readr)
library(treerabid)
library(data.table)
library(dplyr)
library(tidyr)
library(magrittr)
library(foreach)
library(iterators)
library(igraph)
library(doParallel)
source(here("R/boot_trees_simulate_location.R")) 

animal_cases <- read_csv(here("data/clean/dat_animal_for_Romblon_tm_trees.csv")) 
pop_size_by_cell <- read_csv(here("output/population_sizes/Romblon_pop_by_cell.csv"))

case_dates <- data.table(id_case = animal_cases$id,
                         date = animal_cases$date_symptoms_started_proxy)

# Create parameter combinations for different scenarios: using just epi data with & w/out pruning
tree_pars <- tidyr::expand_grid(si_pdist = "lnorm", 
                                dist_pdist = "weibull",
                                convolve = "baseline",
                                prune = TRUE, 
                                time_cutoff = c(1, 0.99),
                                dist_cutoff = c(1, 0.99),
                                use_known = FALSE, 
                                use_gen = c(FALSE, TRUE),
                                location = c("centroid","simulated"),
                                nsim = 1000)

tree_pars %<>%
  # eliminate parameter combination that prunes only by distance
  filter(!(dist_cutoff == 0.99 & time_cutoff == 1)) %>%
  # change "prune" parameter to correctly match the cutoffs
  mutate(prune = ifelse(time_cutoff == 1 & dist_cutoff == 1, FALSE, TRUE))

# Add additional sets for 0.975 and 0.95 pruning cut-offs
tree_pars2 <- tree_pars %>% replace(., tree_pars == 0.99, 0.975) %>% 
  filter(!(dist_cutoff == 1 & time_cutoff == 1)) 
tree_pars3 <- tree_pars %>% replace(., tree_pars == 0.99, 0.95) %>% 
  filter(!(dist_cutoff == 1 & time_cutoff == 1)) 

tree_pars4 <- tidyr::expand_grid(si_pdist = "lnorm", 
                                 dist_pdist = "weibull",
                                 convolve = "baseline",
                                 prune = TRUE, 
                                 time_cutoff = 0.975,
                                 dist_cutoff = 0.99,
                                 use_known = FALSE, 
                                 use_gen = c(FALSE, TRUE),
                                 location = c("centroid","simulated"),
                                 nsim = 1000)

tree_pars <- 
  rbind(tree_pars, tree_pars2, tree_pars3, tree_pars4) %>%
  mutate(scenario = row_number())
rm(tree_pars2, tree_pars3, tree_pars4)
tree_pars$seed <- 23 * 1:nrow(tree_pars) + 8

# Run treerabid to generate and summarise trees, in parallel
cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)

comb <- function(...) {
  mapply('rbind', ..., SIMPLIFY = FALSE, fill = TRUE)
}

trees_all <-
  foreach(i = iter(tree_pars, by = "row"), .combine = comb) %do% {

    if(i$use_gen) {
      lin_dt <- data.table(id_case = animal_cases$id, lineage = animal_cases$lineage)   
    } else {
      lin_dt <- data.table(id_case = animal_cases$id, lineage = 0)      
    }
    
    if(i$location == "centroid"){
      ttrees <- #use normal boot_trees with fixed location
        boot_trees(id_case = animal_cases$id,
                   id_biter = NULL,
                   x_coord = animal_cases$utm_easting_centroid,
                   y_coord = animal_cases$utm_northing_centroid,
                   owned = animal_cases$owned,
                   date_symptoms = animal_cases$date_symptoms_started_proxy,
                   days_uncertain = animal_cases$days_uncertain,
                   use_known_source = FALSE,
                   exclude_progen = animal_cases$exclude_progen, #all FALSE
                   lineages = lin_dt,
                   prune = i$prune,
                   si_fun = treerabid::si_lnorm1,
                   dist_fun = treerabid::dist_weibull1,
                   cutoff = c(time = i$time_cutoff, dist = i$dist_cutoff),
                   params = treerabid::params_treerabid,
                   N = 1000,
                   seed = i$seed)
    } else {
      ttrees <- #use adapted boot_trees with location sampler
        boot_trees_simulate_location(id_case = animal_cases$id,
                                     id_biter = NULL,
                                     locality = animal_cases$locality,
                                     pop_by_cell = pop_size_by_cell,
                                     owned = animal_cases$owned, 
                                     date_symptoms = animal_cases$date_symptoms_started_proxy, 
                                     days_uncertain = animal_cases$days_uncertain, 
                                     exclude_progen = animal_cases$exclude_progen, #all FALSE
                                     lineages = lin_dt,
                                     prune = i$prune,
                                     si_fun = treerabid::si_lnorm1,
                                     dist_fun = treerabid::dist_weibull1,
                                     cutoff = c(time = i$time_cutoff, dist = i$dist_cutoff),
                                     params = treerabid::params_treerabid,
                                     N = 1000,
                                     seed = i$seed)
    }

    # way to join back up with scenarios
    id <- i[, "scenario"]
    
    # Summarize the trees
    links_all <- build_all_links(ttrees, N = i$nsim)
    links_consensus_consistent <- 
      build_consensus_links(links_all, case_dates = case_dates,
                            lineages = lin_dt, known_progens = NULL, 
                            link_all = TRUE)
    
    # consensus without loops fixed
    links_consensus_raw <- links_all[links_all[, .I[which.max(links)], 
                                               by = c("id_case")]$V1]
    
    # majority & mcc
    tree_ids <- c(mcc = 
                    build_consensus_tree(links_consensus_consistent, ttrees, links_all,
                                         type = "mcc", output = "sim"), 
                  majority = 
                    build_consensus_tree(links_consensus_consistent, ttrees, links_all,
                                         type = "majority", output = "sim"))
    ttrees$mcc <- ifelse(ttrees$sim %in% tree_ids["mcc"], 1, 0)
    ttrees$majority <- ifelse(ttrees$sim %in% tree_ids["majority"], 1, 0)
    
    # also get reff
    reff <- ttrees[, .(reff = .N), 
                   by = c("sim", "id_progen")][, 
                                               .(sim, id_case = id_progen, reff)]
    ttrees <- reff[ttrees, on = c("sim", "id_case")]      
    setnafill(ttrees, cols = "reff", fill = 0)
    
    list(ttrees = cbind(ttrees, id), 
         links_consensus_consistent = cbind(links_consensus_consistent, id), 
         links_consensus_raw = cbind(links_consensus_raw, id))
    
  }

parallel::stopCluster(cl)

# Output files ---
fwrite(trees_all$ttrees, here("output/tm_trees/trees_all.gz"))
fwrite(trees_all$links_consensus_consistent, "output/tm_trees/links_consensus_consistent.csv") 
fwrite(trees_all$links_consensus_raw, "output/tm_trees/links_consensus_raw.csv")

# for joining back up with metadata on parameters used (don't want to add cols to 
# the full large data.tables)
fwrite(tree_pars, "output/tm_trees/scenarios.csv")




