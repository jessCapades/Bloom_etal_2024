# input: sizes of ab and p1 cell
# output: eccentricity output for each cell

cell_eccentricity = function(subset_embryo_measurements) {
  
  p1_eccentricity = (sqrt((subset_embryo_measurements$p1_height/2)^2 - (subset_embryo_measurements$p1_length/2)^2)) / (subset_embryo_measurements$p1_height/2)
  ab_eccentricity = (sqrt((subset_embryo_measurements$ab_height/2)^2 - (subset_embryo_measurements$ab_length/2)^2)) / (subset_embryo_measurements$ab_height/2)
  
  cell_roundness   = cbind(p1_eccentricity, ab_eccentricity)
  
  na_positions = which(is.na(cell_roundness))
  na_positions_for_df = na_positions - 45
  na_emb_dimensions = subset_embryo_measurements[na_positions_for_df, c(6,7)]
  na_emb_ab_eccentricity = (sqrt((na_emb_dimensions$ab_length/2)^2 - (na_emb_dimensions$ab_height/2)^2)) / (na_emb_dimensions$ab_length/2)
  
  na_embs_count = 1
  for(i in na_positions){
    
    cell_roundness[i] = na_emb_ab_eccentricity[na_embs_count]
    na_embs_count = na_embs_count +1
  }
  
  
  return(cell_roundness)
}

ab_p1_ratios = function(subset_embryo_measurements) {
  
  p1_ab_perimeter_ratio = subset_embryo_measurements$p1_perimeter / subset_embryo_measurements$ab_perimeter
  p1_ab_area_ratio      = subset_embryo_measurements$p1_area / subset_embryo_measurements$ab_area
  
  cell_dimension_ratios = cbind(p1_ab_perimeter_ratio, p1_ab_area_ratio)
  
  return(cell_dimension_ratios)
}