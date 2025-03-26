#' signal analysis for par signal linescans
#' 
#' normalized x-coordinates out of 1 and creates a rolling average of par signal
#' @param signal_and_coordinates x-coordinates and signal measurements from linescans
#' @param unique_embryo_ids list that contains all unique embryo ids
#' @return returns lists with normalized x_coordinates and y-values with rolling average

signal_calculations = function(signal_and_coordinates, unique_embryo_ids){
  
  avg_par_signal_measurements = list()
  
  for (unique_embryo_id in unique_embryo_ids){
    unique_measurement_timepoint = subset(signal_and_coordinates, signal_and_coordinates$embryo_id == unique_embryo_id)
    rolling_average_par_signal = rollmean(unique_measurement_timepoint$y.value, k = 10, fill= NA)
    avg_par_signal_measurements = append(avg_par_signal_measurements, rolling_average_par_signal)
  }
  signal_and_coordinates$avg_y_value = as.numeric(avg_par_signal_measurements)
  return(signal_and_coordinates)
}


x_coordinates_calculations = function(signal_and_coordinates, unique_embryo_ids){
  
  normalized_par_signal_x_coordinates = list()
  
  for(unique_embryo_id in unique_embryo_ids){
    
    unique_measurement_timepoint = subset(signal_and_coordinates,
                                          signal_and_coordinates$embryo_id == unique_embryo_id)
    x_values_across_embryo = unique_measurement_timepoint$x.value
    maximum_x_value = max(x_values_across_embryo)
    normalized_x_values = x_values_across_embryo/maximum_x_value
    normalized_par_signal_x_coordinates = append(normalized_par_signal_x_coordinates, normalized_x_values)
  }

  signal_and_coordinates$normalized_x_values = as.numeric(normalized_par_signal_x_coordinates)
  return(signal_and_coordinates)
  
}


background_signal_interpolation = function(norm_smoothed_signal_and_coordinates){

  par_background_signal = subset(norm_smoothed_signal_and_coordinates,
                                     norm_smoothed_signal_and_coordinates$signal_or_background == "background")
  
  par_fluorescent_signal = subset(norm_smoothed_signal_and_coordinates,
                                      norm_smoothed_signal_and_coordinates$signal_or_background == "signal")
  
  background_interpolated_x_values = approx(par_background_signal$normalized_x_values, par_background_signal$avg_y_value,
                                            xout = par_fluorescent_signal$normalized_x_values)
  
  par_fluorescent_signal$background_interpolated_values = background_interpolated_x_values$y
  
  par_fluorescent_signal = na.omit(par_fluorescent_signal)
  
  return(par_fluorescent_signal)
}