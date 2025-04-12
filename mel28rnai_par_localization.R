## Visualize par protein localization after knock down of mel-28

library(ggplot2); library(dplyr); library(reshape); library(Hmisc); library(tidyverse);
library(stats); library(EnvStats); library(cowplot); library(zoo); library(utils)

get_file_info = function(trailingOnly = TRUE, asValues = TRUE) {
  data_file_info = commandArgs(TRUE)
  if (length(data_file_info) < 3){
    stop("Provide the file path to functions,
         data file dir, mel28rnai measurements,
         eg. dev/ mel28_measurements_data_dir/ mel-28_linescan_quantifications.ext")
  } else {
    return(data_file_info)
  }
}

mel28rnai_measurement_file_info = get_file_info()

function_directory = mel28rnai_measurement_file_info[1]
early_embryo_measurements_directory = mel28rnai_measurement_file_info[2]
mel28rnai_measurements_filename = mel28rnai_measurement_file_info[3]

source(paste0(function_directory,"select_variables_reformat_long.R"))
source(paste0(function_directory,"create_factors_sample_sizes.R"))
source(paste0(function_directory,"two_cell_size_comparisons.R"))
source(paste0(function_directory,"four_cell_phenotype_frequencies.R"))
source(paste0(function_directory,"select_variables_reformat_long.R"))
source(paste0(function_directory, "par_signal_analysis.R"))

mel28rnai_linescan_measurements = read.csv(paste0(early_embryo_measurements_directory,
                                                  mel28rnai_measurements_filename))
mel28rnai_488_measurements = subset(mel28rnai_linescan_measurements, mel28rnai_linescan_measurements$channel == 488)
mel28rnai_561_measurements = subset(mel28rnai_linescan_measurements, mel28rnai_linescan_measurements$channel == 561)

# create rolling average of par signal
unique_embryo_488_measurements = unique(mel28rnai_488_measurements$embryo_id)
unique_embryo_561_measurements = unique(mel28rnai_561_measurements$embryo_id)

averaged_par_488_signal = signal_calculations(mel28rnai_488_measurements, unique_embryo_488_measurements)
averaged_par_488_signal_and_norm_coordinates = x_coordinates_calculations(averaged_par_488_signal, unique_embryo_488_measurements)

averaged_par_561_signal = signal_calculations(mel28rnai_561_measurements, unique_embryo_561_measurements)
average_par_561_signal_and_norm_coordinates = x_coordinates_calculations(averaged_par_561_signal, unique_embryo_561_measurements)

# subtract background from signal par_488_signal_and_interpolated_background = background_signal_interpolation(averaged_par_488_signal_and_norm_coordinates)
par_488_signal_and_interpolated_background = background_signal_interpolation(averaged_par_488_signal_and_norm_coordinates)
par_561_signal_and_interpolated_background = background_signal_interpolation(average_par_561_signal_and_norm_coordinates)

par_488_signal_and_interpolated_background$normalized_signal = par_488_signal_and_interpolated_background$avg_y_value - par_488_signal_and_interpolated_background$background_interpolated_values
par_561_signal_and_interpolated_background$normalized_signal = par_561_signal_and_interpolated_background$avg_y_value - par_561_signal_and_interpolated_background$background_interpolated_values

par_signal_488_561_combined = bind_rows(par_488_signal_and_interpolated_background, par_561_signal_and_interpolated_background)

pseudocleavage_par_signal = subset(par_signal_488_561_combined, par_signal_488_561_combined$stage == "pseudocleavage")
pn_centration_par_signal = subset(par_signal_488_561_combined, par_signal_488_561_combined$stage == "pn_centration")

# ## create factors for visualization
pseudocleavage_par_signal$embryo_type = factor(pseudocleavage_par_signal$embryo_type, levels = unique(pseudocleavage_par_signal$embryo_type))
pseudocleavage_par_signal$channel = factor(pseudocleavage_par_signal$channel, levels = unique(pseudocleavage_par_signal$channel))
pseudocleavage_par_signal$embryo_id = factor(pseudocleavage_par_signal$embryo_id, levels = unique(pseudocleavage_par_signal$embryo_id))
pn_centration_par_signal$embryo_type = factor(pn_centration_par_signal$embryo_type, levels = unique(pn_centration_par_signal$embryo_type))
pn_centration_par_signal$channel = factor(pn_centration_par_signal$channel, levels = unique(pn_centration_par_signal$channel))
pn_centration_par_signal$embryo_id = factor(pn_centration_par_signal$embryo_id, levels = unique(pn_centration_par_signal$embryo_id))

# ## visualize each stage's par localization
pseudocleavage_vis = ggplot(pseudocleavage_par_signal, mapping = aes(normalized_x_values, normalized_signal, group = embryo_id, color = channel))
pseudocleavage_vis + facet_wrap(~embryo_type + channel) + theme_cowplot(14) + geom_line()
#pseudocleavage_488_vis +  geom_smooth(method = "loess", se = TRUE) + facet_wrap(~embryo_type) + theme_cowplot(14)
ggsave(paste0(early_embryo_measurements_directory,"pseudocleavage_par_localization_mel28rnai_k_20.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


pn_centration_vis = ggplot(pn_centration_par_signal, mapping = aes(normalized_x_values, normalized_signal, group = embryo_id, color = channel))
pn_centration_vis +  geom_line() + facet_wrap(~embryo_type + channel) + theme_cowplot(14)
ggsave(paste0(early_embryo_measurements_directory,"pn_centration_par_localization_mel28rnai_k_20.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)

# think about creating an x range from 0 to 1 with regularly spaced intervals and interpolating all data signal and background to this range
## graph par signal avg and sd

# calculate grouped signal average and sd based on embryo type
# visualize signal average and sd on line graph

