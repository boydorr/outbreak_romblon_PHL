# This function was adapted from a
# [script in the boydorr/PembaRabies GitHub repository](https://github.com/boydorr/PembaRabies/blob/main/6b.trees_main.R)
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

library(ggplot2)
library(data.table)
library(ggforce)
library(animation)
library(dplyr)
library(patchwork)
library(data.table)
library(lubridate)
library(ggthemes)

#' Create series of plots of coloured lineage or transmission chains spreading on a map
#'
#' @param links_consensus Output from `treerabid::build_consensus_links()`.
#' @param case_data A `data.table` with case dates and coordinates (Universal 
#'   Transverse Mercator coordinate system). This needs to include 4 columns, 
#'   named `"id_case"` (with case identifiers), `"symptoms_started"` (with dates 
#'   in `Date` class), `"utm_easting"` (UTM eastings of case locations) and 
#'   `utm_northing"` (UTM northings of case locations).
#' @param map A shape file of the outbreak area, using the Universal Transverse
#'   Mercator coordinate system.
#' @param colour_by Variable by which to colour the cases in the epicurve. One 
#'   of `"lineage_chain"` (default), `"lineage"`, or `"membership"`.
#' @param bezier_frac A number between 0 and 1, indicating the fraction of the
#'   line at which you want the Bezier curve to peak (i.e. 0.5 is halfway).
#' @param bezier_transform Function to create the distance scale for the Bezier 
#'   curve.
#'   
#' @return A series of plots with a coloured epicurve (case histogram) above a 
#'   map, with each plot representing one month and showing the transmissions   
#'   within that month as coloured curves on the map.
animate_trees_on_map <- function(links_consensus, case_data, map,
                                 colour_by = "lineage_chain",
                                 bezier_frac = 0.25,
                                 bezier_transform = function(x) log(1 / x) * 500){
                                   
  if(!(colour_by %in% c("lineage_chain", "lineage", "membership"))){
    stop("`colour_by` must be one of 'lineage_chain', 'lineage', or 'membership'.")
  }
  if(!(is.numeric(bezier_frac) && bezier_frac >= 0 && bezier_frac <= 1)){
    stop("`bezier_frac` must be a number between 0 and 1 (a fraction).")
  }
  
  links_consensus <- 
    links_consensus[case_data[, .(id_case, utm_easting, utm_northing, symptoms_started)], 
                  on = "id_case"] 
  
  links_consensus[[colour_by]] <- factor(links_consensus[[colour_by]])
  
  links_consensus %>%
    select(id_case, x_coord_to = utm_easting, y_coord_to = utm_northing, id_progen, 
           chain_id = {{colour_by}}, date_infectious = symptoms_started) %>%
    mutate(month = floor_date(date_infectious, unit = "month")) %>%
    left_join(select(case_data, id_progen = id_case, x_coord_from = utm_easting,  
                     y_coord_from = utm_northing, date_exposed = symptoms_started)) %>%
    ungroup() %>%
    mutate(group = 1:nrow(.), 
           date_exposed = case_when(is.na(date_exposed) ~ ymd(date_infectious), 
                                    TRUE ~ ymd(date_exposed)), 
           x_coord_from = ifelse(is.na(x_coord_from), x_coord_to, x_coord_from), 
           y_coord_from = ifelse(is.na(y_coord_from), y_coord_to, y_coord_from)) -> chain_segs
  
  
  
  #Get Bezier curves
  bez_pts <- get_bezier_pts(from = as.matrix(select(chain_segs, lat = y_coord_from, long = x_coord_from)),
                            to = as.matrix(select(chain_segs, lat = y_coord_to, long = x_coord_to)),
                            frac = bezier_frac,
                            transform = bezier_transform)
  
  bez_pts <- left_join(chain_segs, bez_pts, by = "group")
  
  bez_pts %>% 
    mutate(date_infectious = floor_date(date_infectious, unit = "month"), 
           date_exposed = floor_date(date_exposed, unit = "month")) -> bez_pts
  # sampled
  bez_pts$sampled <- bez_pts$id_case %in% links_consensus$id_case[links_consensus$lineage != 0]
  
  #Get data for epicurve
  links_consensus %>%
    mutate(month = floor_date(symptoms_started, unit = "months")) %>%
    group_by(month, pick({{colour_by}})) %>%
    summarise(cases = n(),
              start_date = min(month)) -> top_chain_hist
  
  top_chain_hist$chain_id <- top_chain_hist[[colour_by]]
  
  #Colours for plotting
  cols_chains <- c("#CACACA",
                   "#008CF9", "#B80058", "#EBAC23", "#006E00", "#00BBAD",
                   "#D163E6", "#B24502", "#FF9287", "#5954D6", "#00C6F8",
                   "#878500")
  names(cols_chains) <- c(0:9)
  cols_chains <- cols_chains[2:(max(as.numeric(links_consensus[[colour_by]])) + 1)]
  
  #Make gif
  make_gif(bez_pts, top_chain_hist, cols_chains,
           map, 
           date_seqs = seq.Date(from = floor_date(min(links_consensus$symptoms_started),"month"), 
                                to = floor_date(max(links_consensus$symptoms_started),"month"), 
                                by = "month"), 
           fade_out = 1)  
}

# Functions for gif & bezier curves
# see ?ggforce::geom_bezier2 for more details
# from = data.frame with two columns (lat/long)
# to = data.frame with two columns (lat/long)
# frac is where you want the curve to peak along the line (i.e. 0.5 is halfway)
# transform is the function to create the distance scale (you can play with this scaling too)
# adapted from Dudas et al. Curonia
# https://github.com/evogytis/baltic/blob/master/curonia.ipynb
# Output = data.frame with to/mid/from
get_bezier_pts <- function(from, to, frac = 0.8, 
                           transform = function(x) sqrt(1 / x) * 0.5, 
                           min_dist = 50) {
  sign <- ifelse(from[, "long"] > to[, "long"], -1, 1)
  slope <- (to[, "lat"] - from[, "lat"]) / (to[, "long"] - from[, "long"])
  distance <- sqrt((to[, "lat"] - from[, "lat"])^2 + (to[, "long"] - from[, "long"])^2)
  height <- transform(distance)
  
  hdist <- sqrt(height^2 + (distance * frac)^2) # distance between desired height and point along line
  hdist[distance < min_dist] <- 0
  
  from <- data.frame(long = from[, "long"], lat = from[, "lat"], index = 1, group = 1:nrow(from))
  ctrl <- data.frame(
    long = from[, "long"] + hdist * cos(atan(height / distance / frac) + atan(slope)) * sign,
    lat = from[, "lat"] + hdist * sin(atan(height / distance / frac) + atan(slope)) * sign,
    index = 2, group = 1:nrow(from)
  ) # the magic control point
  to <- data.frame(long = to[, "long"], lat = to[, "lat"], index = 3, group = 1:nrow(from))
  
  df <- do.call(rbind, list(from, ctrl, to))
  df %>%
    group_by(group) %>%
    arrange(index) -> df
  return(df) # return coordinate df sorted
}

# Function for making gif
# fade out controls the # of steps 
make_gif <- function(bez_pts, 
                     top_chain_hist,
                     cols_chains, 
                     map, 
                     date_seqs = seq.Date(from = ymd("2002-01-01"), 
                                          to = ymd("2015-01-01"), 
                                          by = "month"), 
                     fade_out = 1) {
  step_alpha <- 0.9/fade_out
  bez_pts$alpha <- 0.9
  infectious <- exposed <- slice(bez_pts, 0)
  
  for (i in seq_len(length(date_seqs))) {
    infectious$alpha <- infectious$alpha - step_alpha
    infectious_now <- filter(bez_pts, date_infectious == date_seqs[i], 
                             index == 1)
    infectious <- bind_rows(infectious_now, 
                            filter(infectious, alpha != 0))
    exposed_now <- filter(bez_pts, id_progen %in% infectious_now$id_case)
    exposed <- bind_rows(filter(exposed, date_infectious >= date_seqs[i]), 
                         exposed_now)
    
    p_map <- ggplot() +
      geom_sf(data = map, fill = "black", color = "black") +
      # lines from infectious to exposed
      geom_bezier2(
        data = exposed_now, 
        aes(x = long, y = lat, group = group, 
            color = chain_id),
        n = 1000,
        alpha = 0.7, 
        size = 2
      ) +
      # exposed case locs
      geom_point(data = filter(exposed, index == 1), 
                 aes(x = x_coord_to, y = y_coord_to, color = chain_id, 
                     shape = is.na(id_progen)), fill = "black", stroke = 1.5) + 
      # infectious case locs
      geom_point(data = infectious, 
                 aes(x = x_coord_to, y = y_coord_to, fill = chain_id, 
                     color = chain_id, 
                     shape = is.na(id_progen), alpha = alpha), size = 3.5) + 
      scale_shape_manual(values = c(21, 22), guide = "none") +
      scale_color_manual(values = cols_chains, guide = "none", 
                         aesthetics = c("color", "fill")) +
      scale_alpha_identity() +
      theme_map() 
    
    p_ts <- 
      ggplot(top_chain_hist) +
      geom_col(aes(x = month, y = cases, fill = chain_id), 
               position = position_stack(reverse = TRUE)) +
      scale_x_date(breaks = date_seqs[i], date_labels = "%B %Y") + 
      scale_fill_manual(values = cols_chains, guide = "none") +
      geom_vline(xintercept = date_seqs[i]) +
      labs(x = "", y = "Number of cases") +
      cowplot::theme_minimal_hgrid() 
    
    p <- p_ts + p_map + plot_layout(heights = c(1, 3))
    print(p)
    ani.pause()
    
  }
}