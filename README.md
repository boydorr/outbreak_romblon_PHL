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
2. Epidemiological data is saved as a csv in the data folder and was analysed with scripts in the R folder, specifically we used:
   - process_outbreak_dat.R for basic epidemiological analysis and description, as well as the files in the R/epi folder
   - for transmission tree inference:
      - process_outbreak_dat_for_Romblon_tm_trees.R to process data
      - run_treerabid.R, the [treerabid](https://github.com/mrajeev08/treerabid/) package (commit [80fb1da](https://github.com/mrajeev08/treerabid/commit/80fb1da8391e764e60975414e17e98e06136a62e)), and an adapted bootstrap function in boot_trees_simulate_location.R to generate transmission trees
      - plot_tm_trees_comparisons.R, plot_tm_trees_epi_gen_sim_loc_pruned.R and helper functions in plot_lineage_ts.R and animate_trees_on_map.R to plot the figures and generate the animation
   - simulate_to_first_case.R for case detection simulations, which uses helper functions defined in simulate_helper_fun.R

 

