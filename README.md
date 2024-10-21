# Combining genomics and epidemiology to investigate a zoonotic outbreak of rabies in Romblon Province, Philippines
Public repository containing all the code and de-identified data for: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13933434.svg)](https://doi.org/10.5281/zenodo.13933434) in *Nature communications* 

Mirava Yuson*, Criselda T Bautista*, Eleanor M Rees, Carlijn Bogaardt, Van Denn D Cruz, Rowan Durrant, Anna Formstone, Daria L Manalo, Duane R Manzanilla, Leilanie Nacion, Hannaniah Aloyon, Jude Karlo Bolivar, Jeromir Bondoc, Christina Cobbold, Mikolaj Kundergorski, Efraim Panganiban, Shynie Vee M Telmo, Jobin Maestro, Mary Elizabeth G Miranda, Nai Rui Chng, Katie Hampson* & Kirstyn Brunker*

*Equal contributions


Analyses were undertaken using R version 4.3.0 (2023-04-21).
Code and analyses take minutes to run. 
Additional details on genetic data and alignments etc undertaken outside of R are detailed in the methods.

## About
Rabies is a viral zoonosis that kills thousands of people annually in low- and middle-income countries across Africa and Asia where domestic dogs are the reservoir. ‘Zero by 30’, the global strategy to end dog-mediated human rabies, promotes a One Health approach underpinned by mass dog vaccination, post-exposure vaccination of bite victims, robust surveillance and community engagement. Using Integrated Bite Case Management (IBCM) and whole genome sequencing (WGS), we enhanced rabies surveillance to detect an outbreak in a formerly rabies-free island province in the Philippines. We inferred that the outbreak was seeded by at least three independent human-mediated introductions that were identified as coming from neighbouring rabies-endemic provinces. Considerable local transmission went undetected, and two human deaths occurred within 6 months of outbreak detection. Suspension of routine dog vaccination due to COVID-19 restrictions likely facilitated rabies spread from these introductions. Emergency response, consisting of awareness measures and ring vaccination, were performed, but swifter and more widespread implementation is needed to contain and eliminate the outbreak and to secure rabies freedom. We conclude that strengthened surveillance making use of new tools such as IBCM, WGS and rapid diagnostic tests can support One Health in action and progress towards the ‘Zero by 30’ goal.

## Methods
This paper combines epidemiological and genomic analyses briefly described below:

### Phylogenetic analysis
1. Trees were generated from publicly available data downloaded from [RABV-GLUE](http://rabv-glue.cvr.gla.ac.uk/#/home) and new sequences generated as part of this study. The sequence data and associated metadata is available in the data folder.
2. Data were processed as follows:
    a) New WGS and Davao WGS (publicly available) combined and aligned:

```bash
cat raw_data/phd_data/ph_wg_ctb_29.fasta raw_data/pgc_data/pgcM_wg_49.fasta > processed_data/data_prep/sequences/ph_concat_wgs_78.fasta 

mafft processed_data/data_prep/sequences/ph_concat_wgs_78.fasta > processed_data/data_prep/sequences/ph_concat_wgs_78.aln.fasta
```

b) Custom R scipts to curate, clean and deduplicate sequences and metadata  
1-curate_sequence_metadata.R  
2_clean_sequence_metadata.R  
3_convert_to_isolateID.R  
4_dedup_seq_data.R  

c) Genome sequences added to existing rabv-glue alignment using mafft -add function

```bash
 mafft --keeplength --add processed_data/data_prep/sequences/ph_concat_wgs_78.aln.fasta processed_data/data_prep/sequences/ph_rabv_glue_wg_isolateIds_dedup.fasta > processed_data/concatenated_alignment/ph_all_combined_690.aln.fasta 
 ```
 
d) Custom R script to merge sequences from the same sample (but submitted as separate GenBank records). Reduced data to 581 sequences.
5_merge_sequences.R  

e) Reconstruct phylogeny using fasttree on 581 sequence dataset

```bash
fasttree -gtr -gamma -nt processed_data/concatenated_alignment/ph_all_merged_by_id_581.fasta > processed_data/trees/ph_all_581_ft.nwk
```

f) Custom R scripts to extract WGS data and root trees (WGS only and all data) by time 
6_wgs_tree.R  

g) Subset tree to WGS only with gotree prune
    
```bash
gotree prune -i processed_data/trees/ph_all_581_ft.nwk -f processed_data/wgs_alignment/wgs.names.txt -r -o processed_data/wgs_alignment/ph_wgs_ft.nwk
```

h) Perform tree dating using R wrapper for lsd2  
7_tree_dating_lsd.R  

i) Use pastml to perform ancestral date reconstruction

Using Philippines administrative level region as state

```bash
pastml -t processed_data/dated_trees/all_wgsrate_lsd_CI.date.nexus -d processed_data/concatenated_alignment/ph_metadata_merged_by_id_581_pastml.csv -s ',' -c Region --prediction_method MPPA --root_date 1909.37 --html_compressed processed_data/pastml_analysis/HTML_compressed_all_mppa_region.html --html processed_data/pastml_analysis/HTML_all_mppa_region.html --upload_to_itol -o processed_data/pastml_analysis/all_mppa_region_pastml --work_dir processed_data/pastml_analysis/all_mppa_region --tip_size_threshold 100
```

Using Philippines administrative level 

```bash
pastml -t processed_data/dated_trees/all_wgsrate_lsd_CI.date.nexus -d processed_data/concatenated_alignment/ph_metadata_merged_by_id_581_pastml.csv -s ',' -c Province --prediction_method MPPA --root_date 1909.37 --html_compressed processed_data/pastml_analysis/HTML_compressed_all_mppa_province.html --html processed_data/pastml_analysis/HTML_all_mppa_province.html --upload_to_itol -o processed_data/pastml_analysis/all_mppa_province_pastml --work_dir processed_data/pastml_analysis/all_mppa_province --tip_size_threshold 100
```


### Transmission tree analysis
1. Epidemiological data is saved as a csv in the data folder and was analysed with scripts in the R folder, specifically we used:
   - process_outbreak_dat.R for basic epidemiological analysis and description, as well as the files in the R/epi folder. Function: Pulls outbreak linelist data from google sheet & does initial data cleaning. Outputs: data/raw/dat_*.csv; data/clean/dat_human.csv; data/clean/dat_animal.csv; data/clean/dat_outbreak_all.csv;
   - process_population_size_data.R has Function: Extracts pop size data and centroid coordinates (UTM) for all 100m2 cells within Romblon province. This Outputs: output/population_sizes/Romblon_pop_by_cell.csv. 
   - patristic_dist.R and process_genetic_data.R. Function: Generation of patristic distance matrix from phylogenetic tree, cluster analysis, and lineage assignment. Outputs: output/phylogenetic_distances/romblonSeq_24_patristicDist_matrix.txt; output/figures/Hamming_distance_frequencies.jpeg; output/phylogenetic_distances/0.002clust.csv. Run both scripts again in case of new genetic data (a new phylogenetic tree); run process_genetic_data.R again in case just want to change the clustering cut-off.

2. for transmission tree inference:
      - process_outbreak_dat_for_Romblon_tm_trees.R to process data. Function: Further data processing and filtering of dat_animal.csv, to prepare data for transmission tree reconstruction. Outputs: data/clean/dat_animal_for_Romblon_tm_trees.csv. Run again if: using a different set of cases (lines 41-45 for date filtering); or changing included columns, dates, lineage assignment, or any other data
      - run_treerabid.R, the [treerabid](https://github.com/mrajeev08/treerabid/) package [v1.0.1](https://zenodo.org/records/13900871), and an adapted bootstrap function in boot_trees_simulate_location.R to generate transmission trees.
      Function: sets up a table with different scenarios (parameter combinations), and runs treerabid functions (boot_trees, build_consensus_links, build_consensus_tree) over these in parallel, to reconstruct transmission trees. Requires boot_trees_simulate_location.R for an adapted boot_trees version, that simulates case locations within each bootstrap. Outputs: output/tm_trees/trees_all.gz; output/tm_trees/links_consensus_raw.csv; output/tm_trees/links_consensus_consistent.csv; output/tm_trees/scenarios.csv. Run again if changes made to the dataset and want to rerun analyses.
      - plot_tm_trees_comparisons.R, plot_tm_trees_epi_gen_sim_loc_pruned.R and helper functions in plot_lineage_ts.R and animate_trees_on_map.R to plot the figures and generate the animation. Function: creates comparative plots with trees for different scenarios. Requires plot_lineage_ts.R for plotting. Outputs: output/figures/consensus_tree_check.jpeg, output/figures/consensus_tree_subset.jpeg. Run again if the underlying transmission trees updated and want plots to reflect latest datasets.
      - plot_tm_trees_epi_gen_sim_loc_pruned.R. Function: creates transmission tree plots [for consensus trees], epicurve, map and animation for a specific scenario (epi + genetic data / simulated locations / pruned by distance and time 0.99). Requires plot_lineage_ts.R for plotting and animate_trees_on_map.R for the animation. Uses simulated case locations from the MCC tree for mapping purposes. Outputs: various plots and animation in output/figures/epi_gen_simulated_locations_prunedDT99. Run again if the underlying transmission trees are updated and want plots to reflect latest datasets/ you want to alter the plots

3. simulate_to_first_case.R for case detection simulations, which uses helper functions defined in simulate_helper_fun.R
   
## License
This repository is licensed under the GPL-2 License. However, it includes code 
adapted from the [treerabid](https://github.com/mrajeev08/treerabid) R package
and the [boydorr/PembaRabies](https://github.com/boydorr/PembaRabies) 
repository, which are licensed under the MIT License, as well as code adapted
from a private repository owned by Kennedy Lushasi, which has been shared with
permission. It also includes a tif file with population density data, cropped
from a [tif by WorldPop](https://hub.worldpop.org/geodata/summary?id=6316), 
which is licensed under CC BY 4.0, and administrative area boundary data
licensed under CC BY 3.0 IGO. For more details, see the `LICENSE` file.




 

