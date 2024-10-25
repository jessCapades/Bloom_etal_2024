# spindle range calculations

calculate_spindle_range = function(spindle_angle_measurements){
  spindle_angle_ranges = list()
  for (embryo in spindle_angle_measurements){
    spindle_range = max(spindle_angle_measurements[3:8,]) - min(spindle_angle_measurements[3:8,])
    spindle_angle_ranges = append(spindle_range)
  }
 
  return(spindle_angle_ranges)
}