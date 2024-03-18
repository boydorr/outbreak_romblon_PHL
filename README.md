# Combining genomics and epidemiology to investigate a zoonotic outbreak: rabies in Romblon Province, Philippines
Public repository for manuscript by: 
Mirava Yuson*, Criselda T Bautista*, Eleanor M Rees, Carlijn Bogaardt, Van Denn D Cruz, Rowan Durrant, Anna Formstone, Daria L Manalo, Duane R Manzanilla, Leilanie Nacion, Hannaniah Aloyon, Jude Karlo Bolivar, Jeromir Bondoc, Christina Cobbold, Mikolaj Kundergorski, Efraim Panganiban, Shynie Vee M Telmo, Jobin Maestro, Mary Elizabeth G Miranda, Nai Rui Chng, Kirstyn Brunker*, Katie Hampson*

*Equal contributions

This repository contains all the code and de-identified data in this study.

Analyses were undertaken using R version 4.3.0 (2023-04-21).
Code and analyses take minutes to run. 
Additional details on genetic data and alignments etc undertaken outside of R are detailed in the methods.

## About
Rabies is a viral zoonotic disease that kills 160 people daily in low/middle-income countries in Africa and Asia where domestic dogs serve as the main vector. ‘Zero by 30’, the global strategy to eliminate human deaths from dog-mediated rabies, promotes a One Health approach underpinned by mass dog vaccination and post-exposure vaccination of bite victims. Using Integrated Bite Case Management (IBCM) and whole genome sequencing (WGS) we enhanced rabies surveillance to detect an outbreak in a formerly rabies-free island province in the Philippines. We inferred that the outbreak was seeded by at least three independent introductions that were identified as coming from nearby rabies-endemic provinces. Considerable transmission went undetected, and two human rabies deaths occurred within 6 months of outbreak detection. We conclude that suspension of routine dog vaccination due to COVID-19 restrictions facilitated rabies spread from these introductions. Emergency response comprising awareness measures and dog vaccination were performed, but swifter and more widespread implementation of these measures are needed to contain and eliminate the outbreak. Strengthened surveillance making use of new tools such as IBCM, WGS  and rapid diagnostic tests can support One Health in action and progress towards the ‘Zero by 30’ goal.

## Methods
This paper combines epidemiological and genomic analyses briefly described below:

1. Trees were generated from publicly available data downloaded from [RABV-GLUE](http://rabv-glue.cvr.gla.ac.uk/#/home) including new sequences generated as part of this study
2. Epidemiological data is saved as a csv in the data folder and was analysed with data in the R folder, specifically we used:
   - process_outbreak_dat.R for basic epidemiological analysis and description, as well as the files in the R/epi folder
   - for transmission tree inference:
      - process_outbreak_dat_for_Romblon_tm_trees.R to process data
      - run_treerabid.R, the [treerabid](https://github.com/mrajeev08/treerabid/) package (commit [80fb1da](https://github.com/mrajeev08/treerabid/commit/80fb1da8391e764e60975414e17e98e06136a62e)), and an adapted bootstrap function in boot_trees_simulate_location.R to generate transmission trees
      - plot_tm_trees_comparisons.R, plot_tm_trees_epi_gen_sim_loc_pruned.R and helper functions in plot_lineage_ts.R and animate_trees_on_map.R to plot the figures and generate the animation
   - simulate_to_first_case.R for case detection simulations, which uses helper functions defined in simulate_helper_fun.R

   

 

