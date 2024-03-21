
# Find closest relatives using a mrca of a set of tips
subset_tree_with_min_additional_tips <- function(tree, mrca, min_additional_tips = 10) {
  # Load required libraries
  library(tidytree)
  
  # Initialize original_num_tips outside the tryCatch block
  original_num_tips <-0
  
  subset_tree <- tryCatch({
    tidytree::tree_subset(tree, mrca, levels_back = 0)
  }, error = function(e) {
    # Set original_num_tips to 1 in case of error
    original_num_tips <<- 1
    
    # If an error occurs, subset the tree at level back 1
    return(tidytree::tree_subset(tree, mrca, levels_back = 1))
  })
  
  # Check if the original number of tips exists and is equal to 1
  if (exists("original_num_tips") && original_num_tips == 1 && length(subset_tree@phylo$tip.label)-original_num_tips >= min_additional_tips) {
      return(subset_tree)
  } else {
    # Count the number of original tips at the selected levels_back
    original_num_tips <- length(subset_tree@phylo$tip.label)

  # Increment levels_back until the tree has at least min_additional_tips on top of the original tips
  levels_back <- 1
  # Start with levels_back = 1
  while (TRUE) {
    # Subset the tree with the incremented levels_back
    subset_tree <- tidytree::tree_subset(tree, mrca, levels_back = levels_back)
    
    # Calculate the number of additional tips
    num_additional_tips <- length(subset_tree@phylo$tip.label) - original_num_tips
    
    # Check if the number of additional tips is at least min_additional_tips
    if (num_additional_tips >= min_additional_tips) {
      break  # Exit the loop if condition is met
    }
    
    # Increment levels_back for the next iteration
    levels_back <- levels_back + 1
    
    # Break if levels_back exceeds the maximum number of levels
    if (levels_back > 10) {
      stop("Levels back exceeded the maximum limit without reaching the desired number of additional tips.")
    }
  }
  
  return(subset_tree)
}
  }

# Example usage:
# subset_tree <- subset_tree_with_min_additional_tips(all_phl.tree, mrca1, min_additional_tips = 10)


plot_reduced_tree <- function(reduced_data, metadata, outbreak_col = "outbreak") {
  # Create the ggtree object
  reduced_plot <- ggtree(reduced_data@phylo, mrsd = '2023-03-01', ladderize = TRUE, size = 0.1, layout = "circular") %<+% metadata
  
  # Add aesthetics and layers
  plot <- reduced_plot +
    aes(color = {{outbreak_col}}) +
    layout_rectangular() +
    scale_color_manual(values = c("black", "darkred", "blue")) +
    guides(colour = "none") +
    theme_tree2() +
    scale_x_continuous(breaks = c(1890, 1920, 1950, 1980, 2020), minor_breaks = seq(1890, 2030, 5)) +
    geom_fruit(geom = geom_tile, mapping = aes(fill = {{outbreak_col}}), width = 5, offset = 0.05) +
    scale_fill_manual(values = c("transparent", "darkred", "blue"),
                      labels = c('', 'Romblon: 2022-23', 'Romblon: 2011-12')) +
    theme(legend.title = element_blank(),
          panel.grid.major = element_line(color = "darkgrey", size = .2, linetype = "dotted"),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
  
  return(plot)
}
