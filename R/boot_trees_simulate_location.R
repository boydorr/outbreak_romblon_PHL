library(doRNG)
library(foreach)

#' Adapted version of treerabid::boot_trees() that samples case locations within
#' given localities, based on population size per 100m2 cell. 
#' N.B. Unlike the original boot_trees, this function assumes use_known_source 
#' is FALSE! (conditional calling of build_known_tree() has been removed)
#'
#' @inheritParams treerabid::build_tree
#' @param N number of trees to build
#' @param seed seed to pass to doRNG to make trees reproducible
#' @param exp_funs functions to export to foreach loop (i.e. for customized
#'  si_fun/dist_fun which relies on an external scripts), may be needed
#'  for certain types of cluster configs
#' @param exp_pkgs packages to export to foreach loop, defaults to
#'  data.table and treerabid, if other dependencies for si_fun or
#'  dist_fun, then pass here. May be needed for
#'  certain types of cluster configs
#'
#' @return a data.table with bootstrapped trees
#' @importFrom foreach foreach
#' @importFrom doRNG %dorng%
boot_trees_simulate_location <- function(id_case,
                       id_biter,
                       locality,
                       pop_by_cell,
                       owned,
                       date_symptoms, # needs to be in a date class
                       days_uncertain,
                       lineages = NULL,
                       exclude_progen = FALSE,
                       #use_known_source = FALSE #Only to be used with use_known_source = FALSE!
                       prune = TRUE,
                       si_fun,
                       dist_fun,
                       cutoff = 0.95,
                       params,
                       min_time = 1e-6,
                       N = 1,
                       seed = 1245,
                       exp_funs = NULL,
                       exp_pkgs = c("data.table", "treerabid", "igraph")) {
  
  if(any(is.na(locality) | is.na(date_symptoms))) {
    stop("Missing data in times or locations!")
  }
  
  # check that date_symptoms is a date class
  if(class(date_symptoms) != "Date") {
    stop("date_symptoms is not of class Date!")
  }
  
  use_known_source <- FALSE
  known_tree <- NULL

  # Check for duplicated id's here
  # (i.e. depending on whether you are using a known source or not)
  msg <- "id_case has duplicated values,
          but you are not using a known set of possible sources for these cases,
          there should only be one record per case id!"
  if(!use_known_source) {
    dups <- tabulate(id_case)
    if(any(dups > 1)) {
      stop(msg)
    }
  } else {
    # any biter ids that are < 1 (i.e. 0, neg numbers, and NAs)
    dups <- tabulate(id_case[id_biter < 1])
    if(any(dups > 1)) {
      stop(msg)
    }
    
  }
  
  foreach(i = seq_len(N),
          .combine = 'rbind', .options.RNG = seed,
          .export = exp_funs,
          .packages = exp_pkgs) %dorng% {
            
            simulated_coordinates <-
              lapply(locality, function(x){
                      idx <- sample(nrow(pop_by_cell[pop_by_cell$loc_ID == x,]), 
                                    size = 1, 
                                    prob = pop_by_cell$population[pop_by_cell$loc_ID == x])
                      pop_by_cell[pop_by_cell$loc_ID == x,
                                  c("utm_easting","utm_northing")][idx,]}) %>%
              purrr::reduce(rbind, make.row.names = FALSE)

            ttree <-
              build_tree(id_case = id_case, id_biter = id_biter,
                         y_coord = simulated_coordinates$utm_northing,
                         x_coord = simulated_coordinates$utm_easting,
                         owned = owned, date_symptoms = date_symptoms,
                         days_uncertain = days_uncertain,
                         exclude_progen = exclude_progen,
                         use_known_source = use_known_source,
                         known_tree = known_tree,
                         lineages = lineages,
                         prune = prune,
                         si_fun,
                         dist_fun,
                         cutoff = cutoff,
                         params = params,
                         min_time = min_time)
            ttree$sim <- i
            ttree
          }
  
}
