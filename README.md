## About the Bloom_etal_2025 repository
This folder consists of code used in the *Caenorhabditis* hybrid incompatibility project that is described here: Hybrid incompatibility emerges at the one-cell stage in interspecies *Caenorhabditis* embryos.  
*Jessica Bloom, Rebecca A. Green, Arshad Desai, Karen Oegema, Scott A. Rifkin*

The code for this project was written in Python 3.1, R 4.3, and Ruby 2.6
***

### Image Editing
Crop embryos and split channels  
*scripts: dic_emb_cropping.ijm, stack_to_slices.ijm, three_channel_emb_cropping.ijm*

Adjust image contrast  
*scripts: split_stacks_to_single_images.ijm, adjust_contrast_for_each_slice.ijm*

### Blinded embryo scoring
Calculate frequency score for each phenotype  
*scripts: blind_embryo_scoring_quantifications_python*

Create a timecourse-ordered heat map from each phenotype frequency score  
*scripts: blinded_evo_phenotypes_scoring.R*

Anonymize embryo image file names  
*scripts: rename-files.rb*

### Early embryo measurements
Visualize measurements made on early embryos during embryogenesis regarding cell shape, size, and spindle dynamics  
*scripts: EHEP_DIC_GLOSv2.R*

### Centrosome Localization and Pronuclear Phenotype Quantification
Visualize frequency of centriole localization phenotypes and quantify pronuclear size at two timepoints in the first cell division  
*scripts: centriole_localization_main.R, centriole_localization_data_handling.R*

### Spindle Phenotype and Centrosome Localization Phenotype Quantification
Visualize frequency of spindle and centriole localization phenotypes  
*scripts: spindle_centrosome_phenotypes_main.R*

### Par Phenotype Frequency Quantification
Visualize frequency of par-2 and par-6 phenotypes after RNAi treatments  
*scripts: par-RNAi_phenotype_quantifications.R*

### Viability and brood size calculations
Calculate brood size and viability on a per worm basis  
*scripts: viability_and_brood_quantifications_python*

### Visualize the daily brood size and viability per species used  
*scripts: viability_brood_counting.R*

### Calculate protein sequence divergence
Calculate the protein sequence divergence scores between species and relative to *C. brenneri*  
Code only written to include species of *Caenorhahbditis*  
*scripts: protein_sequence_divergence_calculations_python*

### Visualize protein sequence divergence
Create a heatmap of clustered protein sequence divergence scores relative to *C. brenneri*   
*scripts: sequence_evolution_heatmap.R*

### Functions used in multiple analyses  
*calculate_spindle_angle_ranges.R, create_factors_sample_sizes.R, four_cell_phenotype_frequencies.R, get_functions.R, select_variables_reformat_long.R, two_cell_size_comparisons.R*

