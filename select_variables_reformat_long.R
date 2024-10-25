
#' subset data
#' 
#' Selects a subset of variables from larger dataset for analysis and reshapes data to long format for 
#' use in ggplot
#' @param all_embryo_metrics dataframe with all data
#' @param embryo_descriptors variable with string data
#' @param embryo_measurements variable with string data
#' @return returns a dataframe with selected data and variables in long format

## split into two functions, one to subset data and one to reshape data


subset_data = function(all_embryo_metrics, embryo_descriptors, embryo_measurements) {
  
  selected_ids_and_measurements = c(embryo_descriptors, embryo_measurements)
  
  subset_embryo_measurements = all_embryo_metrics[, selected_ids_and_measurements]
  
  edited_variable_names = colnames(selected_ids_and_measurements)
  
  subset_embryo_measurements = na.omit(subset_embryo_measurements)
  
  return(subset_embryo_measurements)
}


reshape_data = function(subset_embryo_measurements, embryo_measurements) {
  embryo_data_long = reshape(subset_embryo_measurements,
                                          timevar   = "measurement_type",
                                          times     = embryo_measurements,
                                          v.names   = "measurement_length",
                                          varying   = embryo_measurements,
                                          direction = "long")
  embryo_data_long <- subset(embryo_data_long, embryo_data_long$markers_on != "")
  return(embryo_data_long)
}
