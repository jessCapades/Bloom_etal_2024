#' create_factors_and_sample_sizes
#' 
#' Selects a subset of variables from larger dataset for analysis and reshapes data to long format for 
#' use in ggplot
#' @param reshaped_embryo_measurements long format dataframe that contains the variables used as factors
#' @param data_descriptors list of variables to act as factors from your data frame
#' @return returns a dataframe with selected data and variables in long format


create_factors = function(reshaped_embryo_measurements, data_descriptors){
  
    for (descriptor in data_descriptors){
      reshaped_embryo_measurements[[descriptor]] = 
        factor(reshaped_embryo_measurements[[descriptor]], 
               levels = unique(reshaped_embryo_measurements[[descriptor]]))
    }
}

#' Selects a subset of variables from larger dataset for analysis and reshapes data to long format for 
#' use in ggplot
#' @param subset_embryo_measurements dataframe of selected measurements with no NA
#' @param data_descriptors list of variables to act as factors from your data frame
#' @return returns a dataframe with selected data and variables in long format


calculate_sample_sizes = function(subset_embryo_measurements, measurement){
  sample_sizes = list()
  species_names = subset_embryo_measurements$embryo_type
  unique_species_names = unique(species_names)
  
  for (species_name in unique_species_names){
    species_data = filter(subset_embryo_measurements, 
                          embryo_type == species_name & 
                            is.na(measurement) == FALSE)
    sample_data_number = nrow(species_data)
    sample_sizes = append(sample_sizes, sample_data_number)
  }
  
  names(sample_sizes) = unique_species_names 
  return(sample_sizes)
}

calculate_sample_sizes_with_markers = function(subset_embryo_measurements, measurement){
  sample_sizes = list()
  species_names = subset_embryo_measurements$embryo_type
  unique_species_names = unique(species_names)
  
  for (species_name in unique_species_names){
    species_data = filter(subset_embryo_measurements, 
                          embryo_type == species_name & 
                            is.na(measurement) == FALSE&
                            markers_on != "")
    sample_data_number = nrow(species_data)
    sample_sizes = append(sample_sizes, sample_data_number)
  }
  
  names(sample_sizes) = unique_species_names 
  return(sample_sizes)
}

