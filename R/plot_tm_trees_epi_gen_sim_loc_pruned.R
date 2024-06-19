library(ggplot2)
library(treerabid)
library(lubridate)
library(dplyr)
library(ggraph)
library(igraph)
library(patchwork)
library(sf)
library(ggthemes)
library(ggforce)
library(animation)
library(here)
library(data.table) 
source(here("R/plot_lineage_ts.R")) 
source(here("R/animate_trees_on_map.R")) 

# load in outputs ----
ttrees <- fread(here("output/tm_trees/trees_all.gz"))
links_consensus <- fread(here("output/tm_trees/links_consensus_consistent.csv")) 
scenarios <- fread(here("output/tm_trees/scenarios.csv"))
animal_cases <- fread(here("data/clean/dat_animal_for_Romblon_tm_trees.csv"))
case_dates <- data.table(id_case = animal_cases$id,
                         symptoms_started = animal_cases$date_symptoms_started_proxy)

# filter trees to those for scenario: epi + gen data / simulated locations / pruned w distance and time cut-offs both 0.99
links_consensus <- links_consensus[scenario == 12]
ttrees <- ttrees[scenario == 12]
mcc_tree <- ttrees[mcc == 1]
maj_tree <- ttrees[majority == 1]

# output table with membership, lineage_chain and lineage identifiers for each case ----
links_consensus %>%
  select(id_case, membership, lineage, lineage_chain) %>%
  arrange(id_case) %>%
  fwrite(here("output/figures/epi_gen_simulated_locations_prunedDT99/constree_epigen_prunedDT99_simlocs_chainids.csv"))

# plot transmission trees ----

# plot transmission trees with edges coloured by transmission chain, but 
# with no points for unsampled nodes
# [Like Figure 4B from Pemba paper, but with no sampled cases]
tree_col_by_membership_no_unsampled_nodes <- 
  plot_lin_ts(links_consensus, case_dates,
              colour_by = "lineage_chain",
              colour_edges = TRUE,
              shape_unsampled = NA,
              id_labels = FALSE,
              timeline_breaks = 3)
ggsave(here("output/figures/epi_gen_simulated_locations_prunedDT99/constree_epigen_prunedDT99_simlocs_nonodes.jpeg"), tree_col_by_membership_no_unsampled_nodes, height = 8, width = 8)


# plot transmission trees with edges coloured by transmission chain, with points 
# and id labels for all (incl unsampled) nodes
tree_col_by_membership_id_labels <- 
  plot_lin_ts(links_consensus, case_dates,
              colour_by = "lineage_chain",
              colour_edges = TRUE,
              shape_unsampled = 21,
              id_labels = TRUE,
              timeline_breaks = 3)
ggsave(here("output/figures/epi_gen_simulated_locations_prunedDT99/constree_epigen_prunedDT99_simlocs_idlabels.jpeg"), tree_col_by_membership_id_labels, height = 8, width = 8)


# plot transmission trees with filled circles coloured by transmission chain 
# for unsampled nodes, but with grey edges.
tree_col_by_membership_grey_edges <- 
  plot_lin_ts(links_consensus, case_dates,
              colour_by = "lineage_chain",
              colour_edges = FALSE,
              shape_unsampled = 21,
              id_labels = NA,
              timeline_breaks = 3)
ggsave(here("output/figures/epi_gen_simulated_locations_prunedDT99/constree_epigen_prunedDT99_simlocs_greyedges.jpeg"), tree_col_by_membership_grey_edges, height = 8, width = 8)


# plot epicurve ----

# [Like Figure 4A from Pemba paper]
epicurve <- plot_epicurve(links_consensus, case_dates,
                          colour_by = "lineage_chain",
                          timeline_breaks = 3)
ggsave(here("output/figures/epi_gen_simulated_locations_prunedDT99/epicurve_consensus_epigen_prunedDT99_simlocs.jpeg"), epicurve, height = 8, width = 8)

# plot map with cases ----

# [Like Figure 4C from Pemba paper, but with no sampled cases]

### change headers to match those required by function
case_data <-
  animal_cases %>%
  left_join(mcc_tree[,c("id_case","x_coord","y_coord")], join_by(id == id_case)) %>%
  mutate(utm_easting = x_coord,
         utm_northing = y_coord) 
### transform map coord ref system to UTM & filter shape file to Romblon island
PHL_shp <- st_read(here("data/gis/PHL_village.shp")) %>% st_transform(crs = 32651)
Romblon_shp <- PHL_shp[PHL_shp$NAME_1 == "Romblon",]
Romblon_shp <- 
  Romblon_shp[Romblon_shp$NAME_2 != "Corcuera" & 
                Romblon_shp$NAME_2 != "Banton" & 
                Romblon_shp$NAME_2 != "Concepcion" & 
                Romblon_shp$NAME_2 != "San Fernando" & 
                Romblon_shp$NAME_2 != "Cajidiocan" & 
                Romblon_shp$NAME_2 != "Magdiwang",]
### plot
map_w_cases <- plot_map_w_lin_cases(links_consensus, case_coords = case_data,
                                    map = Romblon_shp,
                                    colour_by = "lineage_chain")
ggsave(here("output/figures/epi_gen_simulated_locations_prunedDT99/map_consensus_epigen_prunedDT99_simlocs.jpeg"), map_w_cases, height = 8, width = 8)


# Reconstruct complete Figure 4 from Pemba paper ----

# Combine epicurve, transmission trees, and map with cases
tm_trees_main_figure <- 
  (epicurve + 
     wrap_elements(tree_col_by_membership_id_labels) +
     map_w_cases) + 
  plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 1, heights = c(1.5, 7, 2.5), guides = "collect") &
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 15))
ggsave(here("output/figures/epi_gen_simulated_locations_prunedDT99/consensus_epicurve_trees_map.jpeg"), tm_trees_main_figure, height = 11, width = 10)


layout <- "
AAAAABBBB
CCCCCBBBB
CCCCCBBBB
"

fig3 <- epicurve + map_w_cases + wrap_elements(tree_col_by_membership_id_labels) + 
  plot_layout(design = layout) + 
  plot_annotation(tag_levels = list("A","B","C")) + 
  theme(plot.tag = element_text(face = "bold", size = 15))
ggsave(here("output/figures/epi_gen_simulated_locations_prunedDT99/fig3.jpeg"), 
       fig3, height = 8, width = 12)

# ALternative
tm_trees_main_figure <-
  ((epicurve / wrap_elements(tree_col_by_membership_id_labels) + plot_layout(heights = c(2,8), widths = 10)) |
     (map_w_cases + plot_layout(heights = 10, widths = 5))) +
  plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 2, widths = c(10,5), guides = "collect") &
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 15))
ggsave(here("output/figures/epi_gen_simulated_locations_prunedDT99/consensus_epicurve_trees_map.jpeg"), tm_trees_main_figure, height = 10, width = 15)

# Clip with epicurve and animated geographic spread of disease over time ----


### change headers to match those required by function
case_data <-
  animal_cases %>%
  mutate(id_case = id, 
         symptoms_started = date_symptoms_started_proxy,
         utm_easting = utm_easting_centroid,
         utm_northing = utm_northing_centroid)

### make gif with animation and save as a movie (more manageable file size and res than gif/pdf)
saveVideo({ 
  animate_trees_on_map(links_consensus, case_data, map = Romblon_shp,
                       colour_by = "lineage_chain",
                       bezier_frac = 0.25,
                       bezier_transform = function(x) log(1 / x) * 500)},
  video.name = here("output/figures/epi_gen_simulated_locations_prunedDT99/tmtree_animation.mp4"),
  img.name = 'trans_plot', title = 'Reconstructed transmission trees', 
  description = 
    c('Reconstructed transmission trees using lognormal SI and weibull dk'),
  interval = 1, nmax = 50, ani.dev = "png", ani.type = "png",
  ani.width = 1500, ani.height = 2000, ani.res = 300)