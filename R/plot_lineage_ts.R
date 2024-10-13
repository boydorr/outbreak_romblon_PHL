# The functions in this file contain code adapted from original scripts in 
# [boydorr/PembaRabies](https://github.com/boydorr/PembaRabies/blob/main/6b.trees_main.R)
# and a private GitHub repository owned by Kennedy Lushasi. The code was adapted
# by Carlijn Bogaardt in January 2024. 
# Sources have been attributed in the function documentation.
# Code from Kennedy Lushasi was shared with permission. Do not redistribute 
# without further consent from the owner.
# Code from the boydorr/PembaRabies repository is licensed under the MIT 
# license, the full text of which is included below:
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

library(treerabid)
library(data.table)
library(ggthemes)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(igraph)
library(ggraph)

#' Plot transmission trees on a timeline
#' 
#' @description 
#' This function was adapted from the original plot_lin_ts() in Kennedy Lushasi's
#' private repository, by Carlijn Bogaardt in January 2024. Used with permission.
#'
#' @param links_consensus Output from `treerabid::build_consensus_links()`.
#' @param case_dates A `data.table` with 2 columns, named `"id_case"` (with case 
#'   identifiers) and `"symptoms_started"` (with dates in `Date` class).
#' @param colour_by Variable by which to colour the transmission chain nodes and
#'   edges. One of `"lineage_chain"` (default), `"lineage"`, or `"membership"`.
#' @param colour_edges Logical indicating whether edges should be coloured by the 
#'   variable selected for `colour_by`. If `FALSE`, edges will be coloured grey.
#' @param shape_unsampled Identifier indicating the shape for unsampled cases in
#'   the plot. The default `NA` results in unsampled cases not being shown as nodes.
#' @param id_labels Logical indicating whether case identifiers should be added
#'   as labels.
#' @param timeline_breaks Positive integer (or numeric coercible to integer) 
#'   indicating the interval in months between labels on the x-axis timeline.
#' @param figure_tag An uppercase letter or other figure tag to be added to the
#'   upper-left corner of the figure.
#'
#' @return A plot with transmission trees, and a timeline as x-axis.
#' 
#' @author Kennedy Lushasi (most code)
#' @author Carlijn Bogaardt (adaptations only)
plot_lin_ts <- function(links_consensus, case_dates, 
                        colour_by = "lineage_chain",
                        colour_edges = TRUE,
                        shape_unsampled = NA,
                        id_labels = FALSE,
                        timeline_breaks = 1,
                        figure_tag = NULL) {
  
  if(!(colour_by %in% c("lineage_chain", "lineage", "membership"))){
    stop("`colour_by` must be one of 'lineage_chain', 'lineage', or 'membership'.")
  }
  if(!(is.numeric(timeline_breaks) && timeline_breaks%%1==0 && timeline_breaks > 0)){
    stop("`timeline_breaks` must be a positive integer-ish number.")
  }
  if(!(is.logical(colour_edges) && length(colour_edges) == 1)){
    stop("`colour_edges` must be a logical of length one (either TRUE or FALSE).")
  }
  if(!(is.logical(id_labels) && length(colour_edges) == 1)){
    stop("`id_labels` must be a logical of length one (either TRUE or FALSE).")
  }
  
  # Lineages over time and space
  links_consensus <- links_consensus[case_dates, on = "id_case"]
  links_consensus[, days_since_start := as.numeric(ymd(symptoms_started) - ymd(min(symptoms_started)))]
  
  links_gr <- get_graph(from = links_consensus$id_progen, to = links_consensus$id_case,
                        attrs = links_consensus[, c("id_case", "lineage", "lineage_chain", "days_since_start")])
  V(links_gr)$membership <- components(links_gr)$membership
  
  grs <- decompose(links_gr)
  
  shps = c(shape_unsampled, rep(22, 9))
  names(shps) <- c(0, 1:9)
  cols_chains <- c("#CACACA",
                   "#008CF9", "#B80058", "#EBAC23", "#006E00", "#00BBAD",
                   "#D163E6", "#B24502", "#FF9287", "#5954D6", "#00C6F8",
                   "#878500")
  names(cols_chains) <- c(0:(length(cols_chains) -1))
  orders <- order(unlist(lapply(grs, function(x) min(V(x)$days_since_start))))
  
  out <-lapply(orders, function(x) {
    gr_now <- grs[[x]]
    
    if(isTRUE(colour_edges)) {
      gr_now <- set_edge_attr(gr_now, colour_by, 
                              value = vertex_attr(gr_now, colour_by)[1])
    } else {
      gr_now <- set_edge_attr(gr_now, colour_by, value = 0)
    }

    out <- ggraph(gr_now, 'dendrogram', height = days_since_start)
    
    if(length(E(gr_now)) > 0){
      out <- out + geom_edge_bend(aes(edge_colour = factor(.data[[colour_by]])), 
                                  alpha = 0.85,
                                  edge_width = 1)}
    
    out <- out+
      geom_node_point(aes(fill = factor(.data[[colour_by]]),
                          shape = factor(lineage)), color = "black", size = 3) +
      coord_flip(clip = "off") +
      scale_edge_color_manual(values = cols_chains, guide = "none") +
      scale_shape_manual(values = shps, guide = "none") +
      scale_fill_manual(values = cols_chains, guide = "none") +
      ylim(c(0, max(links_consensus$days_since_start))) +
      theme(plot.margin = margin(0, 0, 0, 0), panel.background = element_rect(fill = NA))

    if(isTRUE(id_labels)){ 
      out <- out + 
        geom_node_point(aes(col = factor(.data[[colour_by]]),
                            shape = factor(lineage)), fill = "white", stroke = 1.5, size = 7) +
        geom_node_text(aes(label=name), repel=FALSE) +
        scale_color_manual(values = cols_chains, guide = "none")}
    
    return(out) 
  })
  
  axis_temp <- data.frame(x = 0, y = 0)
  date_ticks <- seq.Date(from = ymd(min(links_consensus$symptoms_started)),
                         ymd(max(links_consensus$symptoms_started)), 
                         by = paste(timeline_breaks, "months"))
  date_vals <- as.numeric(date_ticks - ymd(min(links_consensus$symptoms_started)))
  names(date_vals) <- format.Date(floor_date(date_ticks, unit = "months"), "%b %Y")
  
  axis_dt <-
    ggplot(axis_temp) +
    geom_point(aes(x, y), color = NA) +
    scale_x_continuous(breaks = date_vals, labels = names(date_vals),
                       limits = c(0, max(links_consensus$days_since_start))) +
    theme(plot.margin = margin(0, 0, 0, 0),
          panel.background = element_rect(fill = NA),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.line.x = element_line(color = "black"),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  heights <- c(unlist(lapply(grs, function(x) length(V(x))))[orders], 1)
  out <- c(out, list(axis_dt))

  if(!is.null(figure_tag)){
    out[[1]] <- out[[1]] + 
      labs(tag = figure_tag) + 
      theme(plot.tag = element_text(face = "bold"))
  }
  
  plot_ts <- wrap_plots(out, ncol = 1, heights =  heights)
  
  return(plot_ts)
}

#' Plot an epicurve (case histogram) coloured by lineage or transmission chain
#' 
#' @description
#' This function was adapted from code in the 6b.trees_main.R script in the 
#' [boydorr/PembaRabies](https://github.com/boydorr/PembaRabies/blob/main/6b.trees_main.R)
#' GitHub repository by Carlijn Bogaardt in January 2024. The boydorr/PembaRabies
#' repository is licensed under the MIT license.
#' 
#' @param links_consensus Output from `treerabid::build_consensus_links()`.
#' @param case_dates A `data.table` with 2 columns, named `"id_case"` (with case 
#'   identifiers) and `"symptoms_started"` (with dates in `Date` class).
#' @param colour_by Variable by which to colour the cases in the epicurve. One 
#'   of `"lineage_chain"` (default), `"lineage"`, or `"membership"`.
#' @param timeline_breaks Positive integer (or numeric coercible to integer) 
#'   indicating the interval in months between labels on the x-axis timeline.
#' @param figure_tag An uppercase letter or other figure tag to be added to the
#'   upper-left corner of the figure.
#'
#' @return An epicurve (case histogram) plot, with a timeline as x-axis.
#' 
#' @details
#' This function is adapted from the original work in the boydorr/PembaRabies 
#' repository, licensed under the MIT License.
#' 
#' \preformatted{
#' Copyright (c) 2022 Katie Hampson
#' 
#' Permission is hereby granted, free of charge, to any person obtaining a copy
#' of this software and associated documentation files (the "Software"), to deal
#' in the Software without restriction, including without limitation the rights
#' to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
#' of the Software, and to permit persons to whom the Software is furnished to do so,
#' subject to the following conditions:
#' 
#' The above copyright notice and this permission notice shall be included in all
#' copies or substantial portions of the Software.
#'
#' THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#' IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#' FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#' AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#' LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#' OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#' SOFTWARE.
#' }
#' @author Katie Hampson (main code)
#' @author Carlijn Bogaardt (adaptations only)

plot_epicurve <- function(links_consensus, case_dates, 
                          colour_by = "lineage_chain",
                          timeline_breaks = 1,
                          figure_tag = NULL){
  
  if(!(colour_by %in% c("lineage_chain", "lineage", "membership"))){
    stop("`colour_by` must be one of 'lineage_chain', 'lineage', or 'membership'.")
  }
  if(!(is.numeric(timeline_breaks) && timeline_breaks%%1==0 && timeline_breaks > 0)){
    stop("`timeline_breaks` must be a positive integer-ish number.")
  }
  
  cols_chains <- c("#CACACA",
                   "#008CF9", "#B80058", "#EBAC23", "#006E00", "#00BBAD",
                   "#D163E6", "#B24502", "#FF9287", "#5954D6", "#00C6F8",
                   "#878500")
  names(cols_chains) <- c(0:(length(cols_chains) -1))

  links_consensus <- links_consensus[case_dates, on = "id_case"]
  links_consensus[[colour_by]] <- factor(links_consensus[[colour_by]])
  
  links_consensus %>%
    mutate(month = floor_date(symptoms_started, unit = "months")) %>%
    group_by(month, pick({{colour_by}})) %>%
    summarise(cases = n(),
              start_date = min(month)) -> top_chain_hist
  
  breaks <- seq.Date(min(top_chain_hist$month),
                     max(top_chain_hist$month),
                     by = paste(timeline_breaks,"months"))
  
  epicurve <-
    ggplot(top_chain_hist) +
    geom_col(aes(x = month, y = cases, fill = .data[[colour_by]])) +
    scale_x_date(breaks = breaks, 
                 date_labels = "%b %Y") + 
    scale_fill_manual(values = cols_chains, guide = "none") +
    labs(x = "", y = "Number of cases") +
    cowplot::theme_minimal_hgrid()
  
  if(!is.null(figure_tag)){
    epicurve <- epicurve + labs(tag = figure_tag)
  }
  
  return(epicurve)
}

#' Plot a map with cases coloured by lineage or transmission chain
#' 
#' @description
#' This function was adapted from code in the 6b.trees_main.R script in the 
#' [boydorr/PembaRabies](https://github.com/boydorr/PembaRabies/blob/main/6b.trees_main.R)
#' GitHub repository by Carlijn Bogaardt in January 2024. The boydorr/PembaRabies
#' repository is licensed under the MIT license.
#' 
#' @param links_consensus Output from `treerabid::build_consensus_links()`.
#' @param case_coords A `data.table` with case coordinates, using the Universal 
#'   Transverse Mercator coordinate system. This needs to include 3 columns with 
#'   headers `"id"` (case identifiers), `"utm_easting"` (UTM easting of the case
#'   location), `utm_northing"` (UTM northing of the case location).
#' @param colour_by Variable by which to colour the cases on the map. One of
#'   `"lineage_chain"` (default), `"lineage"`, or `"membership"`.
#' @param map A shape file of the outbreak area, using the Universal Transverse
#'   Mercator coordinate system.
#' @param figure_tag An uppercase letter or other figure tag to be added to the
#'   upper-left corner of the figure.
#'
#' @return A map of the outbreak area, with cases coloured by lineage or 
#'   transmission chain.
#'   
#' @details
#' This function is adapted from the original work in the boydorr/PembaRabies 
#' repository, licensed under the MIT License.
#' 
#' \preformatted{
#' Copyright (c) 2022 Katie Hampson
#' 
#' Permission is hereby granted, free of charge, to any person obtaining a copy
#' of this software and associated documentation files (the "Software"), to deal
#' in the Software without restriction, including without limitation the rights
#' to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
#' of the Software, and to permit persons to whom the Software is furnished to do so,
#' subject to the following conditions:
#' 
#' The above copyright notice and this permission notice shall be included in all
#' copies or substantial portions of the Software.
#'
#' THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#' IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#' FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#' AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#' LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#' OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#' SOFTWARE.
#' }
#' @author Katie Hampson (main code)
#' @author Carlijn Bogaardt (adaptations only)
plot_map_w_lin_cases <- function(links_consensus, case_coords,
                                 map,
                                 colour_by = "lineage_chain",
                                 figure_tag = NULL){
  
  if(!(colour_by %in% c("lineage_chain", "lineage", "membership"))){
    stop("`colour_by` must be one of 'lineage_chain', 'lineage', or 'membership'.")
  }

  case_dt <- data.table(case_data)
  links_consensus <- links_consensus[case_coords[, .(id_case = id, utm_easting, utm_northing)], on = "id_case"] 
  
  cols_chains <- c("grey50",
                   "#008CF9", "#B80058", "#EBAC23", "#006E00", "#00BBAD",
                   "#D163E6", "#B24502", "#FF9287", "#5954D6", "#00C6F8",
                   "#878500")
  names(cols_chains) <- c(0:(length(cols_chains) -1))
  cols_chains <- cols_chains[2:(max(as.numeric(links_consensus[[colour_by]])) + 1)]
    
  map_plot <-
    ggplot() +
    geom_sf(data = map, fill = "light grey", color = "light grey") + # Change colour here!
    geom_point(data = links_consensus, #draws coloured filled squares/circles for sampled/unsampled cases
               aes(x = utm_easting, y = utm_northing,
                   color = factor(.data[[colour_by]]),
                   shape = ifelse(lineage == 0, 16, 15),
                   size = ifelse(lineage == 0, 1.5, 2.5),
                   alpha = ifelse(lineage == 0, 0.8, 1)),
               stroke = 1.5) +
    geom_point(data = filter(links_consensus, lineage != 0), #draws black squares around sampled cases
               aes(x = utm_easting, y = utm_northing),
               color = "black", size = 2.5, shape = 0,
               stroke = 1) +
    scale_color_manual(values = cols_chains, labels = names(cols_chains),
                       name = ifelse(colour_by == "membership", "Transmission chain", "Lineage"),
                       drop = TRUE) +
    scale_shape_identity() +
    scale_size_identity() +
    scale_alpha_identity() +
    theme_map() +
    theme(legend.position = "bottom", plot.margin = margin(0, 0, 0, 0))
  
  if(colour_by != "membership"){
    map_plot <- map_plot + guides(color = guide_legend(override.aes = list(shape = 15)))     
  }

  if(!is.null(figure_tag)){
    map_plot <- map_plot + labs(tag = figure_tag) + 
      theme(plot.tag = element_text(face = "bold"))
  }
  
  return(map_plot)
}


