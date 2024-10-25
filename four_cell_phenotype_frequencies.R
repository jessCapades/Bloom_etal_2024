## Four Cell Phenotype Documentation

four_cell_phenotype_frequencies = function(subset_embryo_measurements) {
  subset_embryo_measurements = filter(subset_embryo_measurements, four_cell_division_phenotype != "")
  four_cell_phenotype_categories = unique(subset_embryo_measurements$four_cell_division_phenotype)
  species_types = c("C. elegans", "C. brenneri", "hybrid")
  elegans_four_cell_counts = list()
  brenneri_four_cell_counts = list()
  hybrid_four_cell_counts = list()
  
  species_embryo_counts = list()

  embryo_count = 1
  for (species in species_types) {
    
    embryo_count_list = list()
    
    for (i in four_cell_phenotype_categories) {
      
      phenotype_freq = nrow(filter(subset_embryo_measurements, embryo_type == species  & four_cell_division_phenotype == i))
      embryo_count_list = append(embryo_count_list, phenotype_freq)
    }
    names(embryo_count_list) = four_cell_phenotype_categories
    species_embryo_counts[[embryo_count]] = embryo_count_list
    embryo_count = embryo_count +1
  }
  names(species_embryo_counts) = species_types
  return(species_embryo_counts)
}



calculate_four_cell_sample_size = function(four_cell_division_freq) {
  
  phenotypes = unique(four_cell_division_freq$four_cell_division_phenotype)
  sample_sizes = list()
  
  for (phenotype in phenotypes){
    n_rows = filter(four_cell_division_freq, four_cell_division_phenotype == phenotype)
    sample_N = sum(n_rows$Count)
    sample_sizes = append(sample_sizes, sample_N)  
  }
  names(sample_sizes) = phenotypes
  
  return(sample_sizes)
}