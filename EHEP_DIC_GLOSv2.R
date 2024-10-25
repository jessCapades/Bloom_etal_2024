## ----setup, include=FALSE---------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---------------------------------------------------------------------------------
library(ggplot2); library(dplyr); library(reshape); library(Hmisc); library(tidyverse); library(stats); library(EnvStats); library(cowplot)
## ---------------------------------------------------------------------------------
get_file_info = function(trailingOnly = TRUE, asValues = TRUE) {
  data_file_info = commandArgs(TRUE)
  if (length(data_file_info) < 5){
    stop("Provide the file path to functions,
         data file dir, early emb divsion measurements,
         one cell phenotypes categorizations, and spindle angle measurements ,
         eg. dev/ centriole_data_dir/ centriole_data.ext one_cell_phenotypes.ext spindle_angle_measurements.ext")
  } else {
    return(data_file_info)
  }
}

early_embryogensis_measurement_file_info = get_file_info()
function_directory = early_embryogensis_measurement_file_info[1]
early_embryo_measurements_directory = early_embryogensis_measurement_file_info[2]
early_embryo_measurements_filename = early_embryogensis_measurement_file_info[3]
one_cell_phenotypes_qualifications_filename = early_embryogensis_measurement_file_info[4]
spindle_angle_range_measurements_filename = early_embryogensis_measurement_file_info[5]
## ---------------------------------------------------------------------------------

source(paste0(function_directory,"select_variables_reformat_long.R"))
source(paste0(function_directory,"create_factors_sample_sizes.R"))
source(paste0(function_directory,"two_cell_size_comparisons.R"))
source(paste0(function_directory,"four_cell_phenotype_frequencies.R"))
source(paste0(function_directory,"select_variables_reformat_long.R"))


## ---------------------------------------------------------------------------------
ehep_glos_measurements = read.csv(paste0(early_embryo_measurements_directory,early_embryo_measurements_filename), header = TRUE)


## ---------------------------------------------------------------------------------
eggshell_descriptors = c("embryo_id", "embryo_type", "markers_on")
eggshell_measurements = c("eggshell_length")
eggshell_lengths = subset_data(ehep_glos_measurements, eggshell_descriptors, eggshell_measurements)


## ---------------------------------------------------------------------------------
eggshell_long = reshape_data(eggshell_lengths, eggshell_measurements)


## ---------------------------------------------------------------------------------
data_descriptors = list("embryo_type", "markers_on")
create_factors(eggshell_long, data_descriptors = data_descriptors)

sample_sizes = calculate_sample_sizes(eggshell_lengths,
                                      "eggshell_length")


## ---------------------------------------------------------------------------------
eggshell_length_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                                Hybrid"),
                                "\n(N=",c(sample_sizes$`C. elegans`,
                                        sample_sizes$`C. brenneri`, sample_sizes$hybrid,")",sep=""))

eggshell_length_vis = eggshell_long %>%
                      mutate(embryo_type = 
                               fct_relevel(embryo_type,
                                           "C. elegans",
                                           "C. brenneri",
                                           "hybrid")) %>%
                      ggplot(eggshell_long,
                             mapping = aes(embryo_type,
                                           measurement_length))

eggshell_length_vis + geom_boxplot() + geom_point(position = position_jitter(w= 0.3, h = 0), size = 3, aes(shape = markers_on, group = embryo_type)) + scale_y_continuous(name="Eggshell Length(uM)",limits = c(0, 70)) + stat_n_text() + theme_cowplot(14)


## ---------------------------------------------------------------------------------
ggsave(paste0(early_embryo_measurements_directory,"EHEP_DIC_GLOS_06223_eggshell_lengths.pdf"),width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300 )


## ---------------------------------------------------------------------------------
cell_dist_descriptors = c("embryo_id", "embryo_type", "markers_on")
cell_dist_measurements = c("to_cleavage_length", "between_cell_length")
cell_dist_lengths = subset_data(ehep_glos_measurements, cell_dist_descriptors, cell_dist_measurements)


## ---------------------------------------------------------------------------------
cell_dist_long = reshape_data(cell_dist_lengths, cell_dist_measurements)


## ---------------------------------------------------------------------------------
data_descriptors = list("embryo_type", "markers_on")
create_factors(cell_dist_long, data_descriptors = data_descriptors)

sample_sizes = calculate_sample_sizes(cell_dist_lengths,
                                      "to_cleavage_length")


## ---------------------------------------------------------------------------------
cell_dist_length_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                                Hybrid"),
                                "\n(N=",c(sample_sizes$`C. elegans`,
                                        sample_sizes$`C. brenneri`,                                   sample_sizes$hybrid),")",sep="")

cell_dist_vis = cell_dist_long %>%
                      mutate(embryo_type = 
                               fct_relevel(embryo_type,
                                           "C. elegans",
                                           "C. brenneri",
                                           "hybrid")) %>%
                      ggplot(cell_dist_long,
                             mapping = aes(embryo_type,
                                           measurement_length,
                                           color = measurement_type))

cell_dist_vis + geom_boxplot() + geom_point(position = position_dodge(0.75), aes(shape = markers_on, group = measurement_type), size = 3 ) + scale_y_continuous(name="AB to P1 Separation(uM)",limits = c(0, 35)) + theme(panel.background = element_rect(fill = "white")) + scale_x_discrete(labels= cell_dist_length_Nlab)



## ---------------------------------------------------------------------------------
ggsave(paste0(early_embryo_measurements_directory,"EHEP_DIC_GLOS_06223_cell_dists.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
cell_dist_corr = ggplot(cell_dist_lengths, 
                        mapping=aes(between_cell_length,
                                    to_cleavage_length,
                                    color = embryo_type))

cell_dist_corr + geom_point() + scale_y_continuous(name="AB to P1 Separation(uM)\n To cleavage length",limits = c(0, 35)) + stat_n_text() + theme(panel.background = element_rect(fill = "white"))


## ---------------------------------------------------------------------------------
ggsave(paste0(early_embryo_measurements_directory,"EHEP_DIC_GLOS_06223_cell_dists_corr.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
embryo_descriptors = c("embryo_id", "embryo_type", "markers_on")
spindle_measurements = c("spindle_angle", "spindle_length")
spindle_angles = subset_data(ehep_glos_measurements, embryo_descriptors, spindle_measurements)

spindle_angles$spindle_angle = as.numeric(spindle_angles$spindle_angle) 
spindle_angles = spindle_angles %>% filter(!str_detect(spindle_angle,"^\\s*[0-9]*\\s*$"))


## ---------------------------------------------------------------------------------
spindle_angles$spindle_angle = abs(spindle_angles$spindle_angle)


## ---------------------------------------------------------------------------------
spindle_angles_long = reshape_data(spindle_angles, spindle_measurements)
spindle_angles_long$measurement_length = as.numeric(spindle_angles_long$measurement_length)


## ---------------------------------------------------------------------------------
data_descriptors = list("embryo_type", "markers_on")
create_factors(spindle_angles_long, data_descriptors = data_descriptors)

sample_sizes = calculate_sample_sizes(spindle_angles,
                                      "spindle_angle")



## ---------------------------------------------------------------------------------
spindle_angle_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                                Hybrid"),"\n(N=", c(sample_sizes$'C. elegans',
                                                    sample_sizes$'C. brenneri',
                                                    sample_sizes$hybrid),")",sep="")

spindle_angles_only = filter(spindle_angles_long, measurement_type == "spindle_angle")

spindle_angle_vis = spindle_angles_only %>%
                      mutate(embryo_type = 
                               fct_relevel(embryo_type,
                                           "C. elegans",
                                           "C. brenneri",
                                           "hybrid")) %>%
                ggplot(spindle_angles_only, mapping = aes(embryo_type,
                                           measurement_length))

spindle_angle_vis + geom_boxplot() + geom_point(position = position_jitter(w=0.2, h=0), size = 4, aes(shape = markers_on)) + labs(y = "Spindle Angle from Horizontal(deg)", x = "Cross Type") + scale_x_discrete(labels= spindle_angle_Nlab) + theme_cowplot(14)

## ---------------------------------------------------------------------------------
ggsave(paste0(early_embryo_measurements_directory,"EHEP_DIC_GLOS_06223_spindle_angles.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)

## ---------------------------------------------------------------------------------
embryo_descriptors = c("embryo_id", "embryo_type", "markers_on")
spindle_measurements = c("spindle_length")
spindle_lengths = subset_data(ehep_glos_measurements, embryo_descriptors, spindle_measurements)

spindle_lengths$spindle_length = as.numeric(spindle_lengths$spindle_length) 
spindle_lengths = spindle_lengths %>% filter(!str_detect(spindle_length,"^\\s*[0-9]*\\s*$"))


## ---------------------------------------------------------------------------------
spindle_lengths_long = reshape_data(spindle_lengths, spindle_measurements)


## ---------------------------------------------------------------------------------
data_descriptors = list("embryo_type", "markers_on")
create_factors(spindle_lengths_long, data_descriptors = data_descriptors)

sample_sizes = calculate_sample_sizes(spindle_lengths,
                                      "spindle_length")



## ---------------------------------------------------------------------------------
spindle_lengths_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                            Hybrid"),"\n(N=", c(sample_sizes$`C. elegans`,
                                                sample_sizes$`C. brenneri`,
                                          sample_sizes$hybrid),")",sep="")

spindle_lengths_vis =spindle_lengths_long %>% mutate(embryo_type = 
                               fct_relevel(embryo_type,
                                           "C. elegans",
                                           "C. brenneri",
                                           "hybrid")) %>% ggplot(spindle_lengths_long,
                             mapping = aes(embryo_type,
                                           measurement_length))

spindle_lengths_vis + geom_boxplot() + geom_point(position = position_jitter(w=0.2, h=0), size = 3, alpha = 0.3) + labs(y = "Spindle Length(uM)", x = "Cross Type") + scale_x_discrete(labels= spindle_lengths_Nlab) + scale_y_continuous(name="Spindle Length(uM)",limits = c(0, 30)) + theme(panel.background = element_rect(fill = "white")) + theme_cowplot(14)


## ---------------------------------------------------------------------------------
ggsave(paste0(early_embryo_measurements_directory,"EHEP_DIC_GLOS_06223_spindle_lengths.pdf"), width = 10,
       height = 10, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
embryo_descriptors = c("embryo_id", "embryo_type", "markers_on")
embryo_measurements = c("ab_p_contact_length")
ab_p_contact = subset_data(ehep_glos_measurements, embryo_descriptors = embryo_descriptors, embryo_measurements = embryo_measurements)



## ---------------------------------------------------------------------------------
ab_p_contact_long = reshape_data(ab_p_contact, embryo_measurements)


## ---------------------------------------------------------------------------------
hybrid_ab_p1_lengths = filter(ab_p_contact, embryo_type == "hybrid")
hybrid_ab_p1_media = median(hybrid_ab_p1_lengths$ab_p_contact_length)

cbn_ab_p1_lengths = filter(ab_p_contact, embryo_type == "C. brenneri")
cbn_ab_p1_median = median(cbn_ab_p1_lengths$ab_p_contact_length)

cel_ab_p1_lengths = filter(ab_p_contact, embryo_type == "C. elegans")
cel_ab_p1_median = median(cel_ab_p1_lengths$ab_p_contact_length)


## ---------------------------------------------------------------------------------
data_descriptors = list("embryo_type", "markers_on")
create_factors(ab_p_contact_long, data_descriptors)

sample_sizes = calculate_sample_sizes(ab_p_contact, "ab_p_contact_length")


## ---------------------------------------------------------------------------------
ab_p_contact_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                            Hybrid"),"\n(N=", c(sample_sizes$'C. elegans',
                                                sample_sizes$'C. brenneri',
                                          sample_sizes$hybrid),")",sep="")


ab_p_contact_vis = ab_p_contact_long %>% mutate(embryo_type = 
                               fct_relevel(embryo_type,
                                           "C. elegans",
                                           "C. brenneri",
                                           "hybrid")) %>% ggplot(ab_p_contact_long,
                             mapping = aes(embryo_type,
                                           measurement_length))

ab_p_contact_vis + geom_boxplot() + geom_point(position = position_jitter(w=0.2, h=0), size = 3, alpha = 0.3) + labs(y = "AB:P1 contact length(uM)", x = "Cross Type") + scale_x_discrete(labels= ab_p_contact_Nlab) + scale_y_continuous(limits = c(0, 45)) + theme_cowplot(14)



## ---------------------------------------------------------------------------------
ggsave(paste0(early_embryo_measurements_directory,"EHEP_DIC_GLOS_06223_ab_p_contact.pdf"), width = 10,
       height = 10, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
ab_cell_descriptors = c("embryo_id", "embryo_type", "markers_on")
ab_cell_measurements = c("ab_height", "ab_length")

ab_cell_dimensions = subset_data(ehep_glos_measurements, ab_cell_descriptors, ab_cell_measurements)


## ---------------------------------------------------------------------------------
ab_cell_long = reshape_data(ab_cell_dimensions, ab_cell_measurements)


## ---------------------------------------------------------------------------------
create_factors(ab_cell_long, ab_cell_descriptors)

## not the same number of ab as p1 cells
sample_sizes = calculate_sample_sizes(ab_cell_dimensions, "ab_cell_height")


## ---------------------------------------------------------------------------------
ab_cell_dimensions_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                            Hybrid"),"\n(N=", c(sample_sizes$'C. elegans',
                                                sample_sizes$'C. brenneri',
                                          sample_sizes$hybrid),")",sep="")


ab_cell_dimensions_vis = ab_cell_long %>%  mutate(embryo_type = 
                               fct_relevel(embryo_type,
                                           "C. elegans",
                                           "C. brenneri",
                                           "hybrid")) %>% ggplot(ab_cell_long,
                             mapping = aes(embryo_type,
                                           measurement_length))

ab_cell_dimensions_vis + geom_boxplot() + geom_point(position = position_jitter(w=0.3, h=0), size = 3, aes(shape = markers_on)) + labs(y = "AB:P1 contact length(uM)", x = "Cross Type") + facet_wrap(~measurement_type) + scale_x_discrete(labels= ab_cell_dimensions_Nlab) + scale_y_continuous(limits = c(0, 45)) + theme_cowplot(14)


## ---------------------------------------------------------------------------------
ggsave(paste0(early_embryo_measurements_directory,"EHEP_DIC_GLOS_06223_ab_cell_dim.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
p1_cell_descriptors = c("embryo_id", "embryo_type", "markers_on")
p1_cell_measurements = c("p1_height", "p1_length")

p1_cell_dimensions = subset_data(ehep_glos_measurements, p1_cell_descriptors, p1_cell_measurements)


## ---------------------------------------------------------------------------------
p1_cell_long = reshape_data(p1_cell_dimensions, p1_cell_measurements)


## ---------------------------------------------------------------------------------
create_factors(p1_cell_long, p1_cell_descriptors[-1])

sample_sizes = calculate_sample_sizes(p1_cell_dimensions, "p1_height")


## ---------------------------------------------------------------------------------
p1_cell_dimN = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                            Hybrid"),"\n(N=", c(sample_sizes$`C. elegans`,
                                                sample_sizes$`C. brenneri`,
                                          sample_sizes$hybrid),")",sep="")


p1_cell_dim_vis = p1_cell_long %>%  mutate(embryo_type = 
                               fct_relevel(embryo_type,
                                           "C. elegans",
                                           "C. brenneri",
                                           "hybrid")) %>% ggplot(p1_cell_long,
                             mapping = aes(embryo_type,
                                           measurement_length))

p1_cell_dim_vis + geom_boxplot() + geom_point(position = position_jitter(w=0.2, h=0), size = 3, aes(shape = markers_on)) + labs(y = "Lenght(uM)", x = "Cross Type") + scale_x_discrete(labels= p1_cell_dimN) + scale_y_continuous(limits = c(0, 45)) + theme_cowplot(14) + facet_wrap(~measurement_type)


## ---------------------------------------------------------------------------------
ggsave(paste0(early_embryo_measurements_directory,"EHEP_DIC_GLOS_06223_p1_cell_dim.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
cell_size_ratios_descriptors = c("embryo_id", "embryo_type", "markers_on")
cell_size_measures = c("p1_height", "p1_length", "ab_height", "ab_length", "p1_area", "p1_perimeter", "ab_area", "ab_perimeter")

cell_sizes = subset_data(ehep_glos_measurements, cell_size_ratios_descriptors, cell_size_measures)


## ---------------------------------------------------------------------------------
ab_p1_eccentricity = cell_eccentricity(cell_sizes)
cell_eccentricity = cbind(cell_sizes[cell_size_ratios_descriptors], ab_p1_eccentricity)

p1_ab_ratios = ab_p1_ratios(cell_sizes)
cell_dim_ratios = cbind(cell_sizes[cell_size_ratios_descriptors], p1_ab_ratios)


## ---------------------------------------------------------------------------------
hybrid_area_ratios = filter(cell_dim_ratios, embryo_type == "hybrid")
hybrid_area_ratios_median = median(hybrid_area_ratios$p1_ab_area_ratio)

cbn_area_ratios = filter(cell_dim_ratios, embryo_type == "C. brenneri")
cbn_area_ratios_median = median(cbn_area_ratios$p1_ab_area_ratio)

cel_area_ratios = filter(cell_dim_ratios, embryo_type == "C. elegans")
cel_area_ratios_median = median(cel_area_ratios$p1_ab_area_ratio)


## ---------------------------------------------------------------------------------
cell_sizes_long = reshape_data(cell_sizes, cell_size_measures[-c(1:4)])

cell_eccentricity_long = reshape_data(cell_eccentricity, colnames(ab_p1_eccentricity))

p1_ab_ratios_long = reshape_data(cell_dim_ratios, colnames(p1_ab_ratios))


## ---------------------------------------------------------------------------------
cell_dim_descriptors = c("embryo_type", "markers_on")
cell_dim_measures = list(cell_sizes_long, cell_eccentricity_long, p1_ab_ratios_long)

for (cell_dim_measure in cell_dim_measures){
  
  create_factors(cell_dim_measure, cell_dim_descriptors)
}

cell_size_N = calculate_sample_sizes(cell_sizes, "p1_height")
cell_eccentricity_N = calculate_sample_sizes(cell_eccentricity, "p1_eccentricity")
cell_dim_ratio_N = calculate_sample_sizes(cell_dim_ratios, "p1_ab_perimeter_ratio")


## ---------------------------------------------------------------------------------
cell_size_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                            Hybrid"),"\n(N=", c(cell_size_N$'C. elegans',
                                                cell_size_N$'C. brenneri',
                                          cell_size_N$hybrid),")",sep="")


cell_perimeter_vis = cell_sizes_long[grep("perimeter", cell_sizes_long$measurement_type), ] %>% mutate(embryo_type = 
                               fct_relevel(embryo_type,
                                           "C. elegans",
                                           "C. brenneri",
                                           "hybrid")) %>% ggplot(cell_sizes_long[grep("perimeter", cell_sizes_long$measurement_type), ],
                                                      mapping = aes(embryo_type, 
                                                      measurement_length))

cell_perimeter_vis + geom_boxplot() + geom_point(position = position_jitter(w=0.3, h=0), size = 3, aes(shape = markers_on)) + labs(y = "Perimeter(uM)", x = "Cross Type") + scale_x_discrete(labels= cell_size_Nlab) + scale_y_continuous(limits = c(0, 150)) + theme_cowplot(14) + facet_wrap(~measurement_type)


## ---------------------------------------------------------------------------------
ggsave(paste0(early_embryo_measurements_directory,"EHEP_DIC_GLOS_06223_two_cell_perimeter.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
cell_size_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                            Hybrid"),"\n(N=", c(cell_size_N$`C. elegans`,
                                                cell_size_N$`C. brenneri`,
                                          cell_size_N$hybrid),")",sep="")


cell_area_vis = cell_sizes_long[grep("area", cell_sizes_long$measurement_type), ] %>% mutate(embryo_type = 
                               fct_relevel(embryo_type,
                                           "C. elegans",
                                           "C. brenneri",
                                           "hybrid")) %>% ggplot(cell_sizes_long[grep("area", cell_sizes_long$measurement_type), ],
                                                      mapping = aes(embryo_type, 
                                                      measurement_length))

cell_area_vis + geom_boxplot() + geom_point(position = position_jitter(w=0.2, h=0), size = 3, aes(shape = markers_on)) + labs(y = expression("Area"(uM^2)), x = "Cross Type") + scale_x_discrete(labels= cell_size_Nlab) + scale_y_continuous(limits = c(0, 1300)) + theme_cowplot(14) + facet_wrap(~measurement_type)


## ---------------------------------------------------------------------------------
ggsave(paste0(early_embryo_measurements_directory,"EHEP_DIC_GLOS_06223_two_cell_area.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
cell_eccentricity_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                            Hybrid"),"\n(N=", c(cell_eccentricity_N$'C. elegans',
                                                cell_eccentricity_N$'C. brenneri',
                                          cell_eccentricity_N$hybrid),")",sep="")


cell_eccentricity_vis = cell_eccentricity_long %>% mutate(embryo_type = 
                               fct_relevel(embryo_type,
                                           "C. elegans",
                                           "C. brenneri",
                                           "hybrid")) %>% ggplot(cell_eccentricity_long,
                                mapping = aes(embryo_type, 
                                measurement_length))

cell_eccentricity_vis + geom_boxplot() + geom_point(position = position_jitter(w=0.2, h=0), size = 3, aes(shape = markers_on)) + labs(y = "Eccentricity", x = "Cross Type") + scale_x_discrete(labels= cell_eccentricity_Nlab) + scale_y_continuous(limits = c(0, 1)) + theme_cowplot(14) + facet_wrap(~measurement_type)

## ---------------------------------------------------------------------------------
ggsave(paste0(early_embryo_measurements_directory,"EHEP_DIC_GLOS_06223/EHEP_DIC_GLOS_06223_eccentricity.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
cell_ratio_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                            Hybrid"),"\n(N=", c(cell_dim_ratio_N$'C. elegans',
                                                cell_dim_ratio_N$'C. brenneri',
                                          cell_dim_ratio_N$hybrid),")",sep="")


cell_dim_ratio_vis = p1_ab_ratios_long %>%
                      mutate(embryo_type = 
                               fct_relevel(embryo_type,
                                           "C. elegans",
                                           "C. brenneri",
                                           "hybrid")) %>%ggplot(p1_ab_ratios_long[grep("area", p1_ab_ratios_long$measurement_type), ],
                                mapping = aes(embryo_type, 
                                measurement_length))

cell_dim_ratio_vis + geom_boxplot() + geom_point(position = position_jitter(w=0.3, h=0), size = 2, alpha = 0.3) + labs(y = "P1:AB area", x = "Cross Type") + scale_x_discrete(labels= cell_ratio_Nlab) + scale_y_continuous(limits = c(0, 1.3)) + theme_cowplot(14) + facet_wrap(~measurement_type) + scale_color_grey()



## ---------------------------------------------------------------------------------

ggsave(paste0(early_embryo_measurements_directory,"EHEP_DIC_GLOS_06223_ab_perim_ratio.pdf"), width = 10,
       height = 10, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
p2_stacking_descriptors = c("embryo_id", "embryo_type", "markers_on")
p2_stacking_measurements = c("p2_stacking")

p2_stacking_dist = subset_data(ehep_glos_measurements, p2_stacking_descriptors, p2_stacking_measurements)



## ---------------------------------------------------------------------------------
p2_stacking_dist_long = reshape_data(p2_stacking_dist, p2_stacking_measurements)


## ---------------------------------------------------------------------------------
data_descriptors = list("embryo_type", "markers_on")
create_factors(p2_stacking_dist_long, data_descriptors)

sample_sizes = calculate_sample_sizes(p2_stacking_dist, "p2_stacking")


## ---------------------------------------------------------------------------------
p2_stacking_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                            Hybrid"),"\n(N=", c(sample_sizes$'C. elegans',
                                                sample_sizes$'C. brenneri',
                                          sample_sizes$hybrid),")",sep="")


p2_stacking_vis = p2_stacking_dist_long %>%
                      mutate(embryo_type = 
                               fct_relevel(embryo_type,
                                           "C. elegans",
                                           "C. brenneri",
                                           "hybrid")) %>%ggplot(p2_stacking_dist_long,
                                mapping = aes(embryo_type, 
                                measurement_length))

p2_stacking_vis + geom_boxplot() + geom_point(position = position_jitter(w=0.2, h=0), size = 3, aes(shape = markers_on)) + labs(y = "P2 to ABp Stacking Distance(uM)", x = "Cross Type") + scale_x_discrete(labels= p2_stacking_Nlab) + scale_y_continuous(limits = c(0, 35)) + theme_cowplot(14) + facet_wrap(~measurement_type)
 


## ---------------------------------------------------------------------------------
ggsave(paste0(early_embryo_measurements_directory,"EHEP_DIC_GLOS_06223_p2_stacking.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
p2_contact_descriptors = c("embryo_id", "embryo_type", "markers_on")
p2_contact_measurements = c("p2_abp_contact_length", "p2_ems_contact_length")

p2_contact_lengths = subset_data(ehep_glos_measurements, p2_contact_descriptors, p2_contact_measurements)


## ---------------------------------------------------------------------------------
p2_contact_lengths_long = reshape_data(p2_contact_lengths, p2_contact_measurements)


## ---------------------------------------------------------------------------------
data_descriptors = list("embryo_type", "markers_on")
create_factors(p2_contact_lengths_long, data_descriptors)

sample_sizes = calculate_sample_sizes(p2_contact_lengths, "p2_abp_contact_length")


## ---------------------------------------------------------------------------------
p2_contact_lengths_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                            Hybrid"),"\n(N=", c(sample_sizes$'C. elegans',
                                                sample_sizes$'C. brenneri',
                                          sample_sizes$hybrid),")",sep="")


p2_contact_vis = p2_contact_lengths_long %>%
                      mutate(embryo_type = 
                               fct_relevel(embryo_type,
                                           "C. elegans",
                                           "C. brenneri",
                                           "hybrid")) %>%ggplot(p2_contact_lengths_long,
                                mapping = aes(embryo_type, 
                                measurement_length))

p2_contact_vis + geom_boxplot() + geom_point(position = position_jitter(w=0.2, h=0), size = 3, aes(shape = markers_on)) + labs(y = "P2 contact Lengths(uM)", x = "Cross Type") + scale_x_discrete(labels= p2_stacking_Nlab) + scale_y_continuous(limits = c(0, 35)) + theme_cowplot(14) + facet_wrap(~measurement_type)



## ---------------------------------------------------------------------------------
ggsave(paste0(early_embryo_measurements_directory,"EHEP_DIC_GLOS_06223_p2_contact_lengths.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
four_cell_descriptors = c("embryo_id", "embryo_type", "markers_on")
four_cell_measurements = c("four_cell_division_phenotype")

four_cell_phenotypes = subset_data(ehep_glos_measurements, four_cell_descriptors, four_cell_measurements)


## ---------------------------------------------------------------------------------
four_cell_phenotype_freq = four_cell_phenotype_frequencies(four_cell_phenotypes)
four_cell_phenotype_frequencies = as.data.frame(do.call(rbind, four_cell_phenotype_freq))


## ---------------------------------------------------------------------------------

four_cell_phenotype_group = four_cell_phenotypes %>%
group_by(embryo_type, four_cell_division_phenotype) %>%
  summarise(Count = n())

four_cell_phenotype_group = filter(four_cell_phenotype_group, four_cell_division_phenotype != "")


## ---------------------------------------------------------------------------------
four_cell_phenotype_group$embryo_type = factor(four_cell_phenotype_group$embryo_type, levels=c("C. elegans", "C. brenneri", "hybrid"))

four_cell_sample_sizes = calculate_sample_sizes(four_cell_phenotypes, "four_cell_division_phenotypes")


## ---------------------------------------------------------------------------------
four_cell_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                            Hybrid"),"\n(N=", c(four_cell_sample_sizes$'C. elegans',
                                                four_cell_sample_sizes$'C. brenneri',
                                          four_cell_sample_sizes$hybrid),")",sep="")


## ---------------------------------------------------------------------------------
embryo_4_cell_division_vis = ggplot(four_cell_phenotype_group, mapping = aes(fill = four_cell_division_phenotype, x = embryo_type, y = Count))



## ---------------------------------------------------------------------------------
embryo_4_cell_division_vis + geom_bar(position="fill", stat="identity") + xlab("Cross Type") + ylab("Fraction of Total Embryos") + theme(text = element_text(size = 15)) + theme_cowplot(16) +scale_x_discrete(labels= four_cell_Nlab)

## ---------------------------------------------------------------------------------
ggsave(paste0(early_embryo_measurements_directory,"EHEP_DIC_GLOS_06223_four_cell_phenotypes.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
four_cell_phenotypes_regrouped = four_cell_phenotypes
 
four_cell_phenotypes_regrouped_summarized = four_cell_phenotypes_regrouped %>%
  group_by(embryo_type, four_cell_division_phenotype, markers_on)%>%
  summarise(Count = n())
four_cell_phenotypes_total_counts = four_cell_phenotypes_regrouped %>%
  group_by(embryo_type, four_cell_division_phenotype)%>%
  summarise(all_embs_count=n())


four_cell_phenotypes_regrouped_summarized = filter(four_cell_phenotypes_regrouped_summarized, four_cell_division_phenotype != "" & markers_on != "")
four_cell_phenotypes_total_counts = filter(four_cell_phenotypes_total_counts, four_cell_division_phenotype != "")

four_cell_phenotype_markers_on_hybrids = filter(four_cell_phenotypes_total_counts, embryo_type == "hybrid")



## ---------------------------------------------------------------------------------
four_cell_phenotypes_with_markers = subset(four_cell_phenotypes_regrouped_summarized, four_cell_phenotypes_regrouped_summarized$embryo_type != "C. brenneri")


## ---------------------------------------------------------------------------------
markers_four_cell_N = calculate_four_cell_sample_size(four_cell_phenotypes_with_markers)

four_cell_phenotypes_with_markers$markers_on <- factor(four_cell_phenotypes_with_markers$markers_on, levels=c("Y", "N"))


## ---------------------------------------------------------------------------------
markers_on_four_cell_Nlab = paste(c("C. el wt", "C. el other", "wild-type", "other", "par-2", "par-6"),
                              "\n(N=", c(four_cell_phenotype_frequencies[1,3],
                                         four_cell_phenotype_frequencies[1,1],
                                         four_cell_phenotype_frequencies[3,3],
                                         four_cell_phenotype_frequencies[3,1],
                                         four_cell_phenotype_frequencies[3,4],
                                         four_cell_phenotype_frequencies[3,2]),")",sep="")


## ---------------------------------------------------------------------------------
markers_on_cell_division_vis = ggplot(four_cell_phenotypes_with_markers, mapping = aes(x = four_cell_division_phenotype, y = Count, fill = markers_on))



## ---------------------------------------------------------------------------------
markers_on_cell_division_vis + geom_bar(position= "dodge", stat="identity")+ xlab("four cell outcome") + ylab("Fraction of Total Embryos")  + theme_cowplot(16) + facet_wrap(~embryo_type) #+ scale_x_discrete(labels= hybrid_four_cell_Nlab)



## ---------------------------------------------------------------------------------
ggsave(paste0(early_embryo_measurements_directory,"EHEP_DIC_GLOS_06223_four_cell_markers_on_phenotypes.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
hybrid_four_cell_categories = subset(four_cell_phenotypes, embryo_type == "hybrid" & markers_on != "" & four_cell_division_phenotype != "")


## ---------------------------------------------------------------------------------
hybrid_four_cell_scatter = ggplot(hybrid_four_cell_categories, aes(markers_on, four_cell_division_phenotype))

hybrid_four_cell_scatter + geom_point(position = position_jitter(w=0.3, h=0), size = 3, alpha = 0.6, stroke = NA) + theme_cowplot(14) + theme(panel.grid.major.y = element_line(color = "red",
                                        size = 0.25,
                                        linetype = 2))


## ---------------------------------------------------------------------------------
ggsave(paste0(early_embryo_measurements_directory,"EHEP_DIC_GLOS_06223_four_cell_hybrid_markers_grid_no_height_jitter_ygrid_line.pdf"), width = 10,
       height = 10, units = "cm",
       device = "pdf", dpi = 300)



## ---------------------------------------------------------------------------------
one_cell_phenotypes = read.csv(one_cell_phenotypes_qualifications_filename, header = TRUE)

## ---------------------------------------------------------------------------------

source(paste0(function_directory,"select_variables_reformat_long.R"))
source(paste0(function_directory,"create_factors_sample_sizes.R"))
source(paste0(function_directory,"two_cell_size_comparisons.R"))
source(paste0(function_directory,"four_cell_phenotype_frequencies.R"))
source(paste0(function_directory,"select_variables_reformat_long.R"))


## ---------------------------------------------------------------------------------
nuclear_phenotypes = one_cell_phenotypes[,c(1:2,3)]

nuclear_phenotypes = subset(nuclear_phenotypes, nuclei_normal.y.n. != "",
                            select = c(embryo_id, embryo_type, nuclei_normal.y.n.))

nuclear_phenotypes_grouped = nuclear_phenotypes %>%
  group_by(embryo_type, nuclei_normal.y.n.)%>%
  summarise(Count = n())

embryo_descriptors = list("embryo_type","nuclei_normal.y.n.")

create_factors(nuclear_phenotypes_grouped, embryo_descriptors)

nuclear_phenotype_sample_size = calculate_sample_sizes(nuclear_phenotypes, "nuclei_normal.y.n.")

nuclear_phenotypes_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                            Hybrid"),"\n(N=", c(nuclear_phenotype_sample_size$jk574,
                                                nuclear_phenotype_sample_size$lkc28,
                                          nuclear_phenotype_sample_size$hybrid),")",sep="")

nuclear_phenotypes_visualization = ggplot(nuclear_phenotypes_grouped,
                                          mapping = aes(x = embryo_type, y = Count, fill = nuclei_normal.y.n.))

nuclear_phenotypes_visualization + geom_bar(position="fill", stat="identity") + xlab("Cross Type") + ylab("Fraction") + theme_cowplot(16) + scale_x_discrete(labels = nuclear_phenotypes_Nlab)

ggsave(paste0(early_embryo_measurements_directory,"EHEP-DIC-GLOS-one_cell_nuclear_phenotypes.pdf"), width = 20, height = 20, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
pn_sizes = one_cell_phenotypes[,c(1,2,4:8)]

pn_sizes = drop_na(pn_sizes)

pn_1_areas = pi*(pn_sizes$pn1_height/2)*(pn_sizes$pn.1_length/2)
pn_2_areas = pi*(pn_sizes$pn2_height/2)*(pn_sizes$pn2_length/2)

pn_1_pn_2_ratio = pn_1_areas/pn_2_areas

pn_sizes$pn_1_areas = pn_1_areas
pn_sizes$pn_2_areas = pn_2_areas
pn_sizes$pn_1_pn_2_ratio = pn_1_pn_2_ratio



## ---------------------------------------------------------------------------------
pn_sample_size = calculate_sample_sizes(pn_sizes, "pn1_height")

pn_data_descriptors = list("embryo_type")
create_factors(pn_sizes, pn_data_descriptors)


## ---------------------------------------------------------------------------------
pn_areas = pn_sizes[, c(1,2,7:10)]

pn_areas_long = reshape(pn_areas, 
        direction = "long",
        varying = c("pn_1_areas","pn_2_areas","pn_1_pn_2_ratio"),
        v.names = "area_measurement",
        timevar = "measurement_type",
        times = c("pn_1_areas","pn_2_areas","pn_1_pn_2_ratio"))

pn_area_data_descriptors = list("embryo_type", "measurement_type")
create_factors(pn_areas_long, pn_data_descriptors)


## ---------------------------------------------------------------------------------
pn_area_ratios = subset(pn_areas_long, measurement_type == "pn_1_pn_2_ratio")

pn_areas_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                            Hybrid"),"\n(N=", c(pn_sample_size$jk574,
                                                pn_sample_size$lkc28,
                                          pn_sample_size$hybrid),")",sep="")

pn_area_ratio_vis = pn_area_ratios %>% mutate(embryo_type = 
                               fct_relevel(embryo_type,
                                           "jk574",
                                           "lkc28",
                                           "hybrid")) %>%ggplot(pn_area_ratios, mapping = aes(embryo_type, area_measurement))

pn_area_ratio_vis + geom_boxplot() + geom_point(position = position_jitter(w=0.2, h=0), size = 5, alpha = 0.6) + theme_cowplot(14) + scale_x_discrete(labels= pn_areas_Nlab) + scale_y_continuous(limits = c(0,2))

## ---------------------------------------------------------------------------------
cbn_pn_ratio_med = median(subset(pn_area_ratios, embryo_type == "lkc28")$area_measurement)
cel_pn_ratio_med = median(subset(pn_area_ratios, embryo_type == "jk574")$area_measurement)
hyb_pn_ratio_med = median(subset(pn_area_ratios, embryo_type == "hybrid")$area_measurement)


## ---------------------------------------------------------------------------------
ggsave(paste0(early_embryo_measurements_directory,"EHEP-DIC-GLOS-pn_area_ratios.pdf"), width = 20, height = 20, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
cytokinesis_phenotypes = one_cell_phenotypes[,c(1,2,9)]
cytokinesis_phenotypes = subset(cytokinesis_phenotypes, cytokinesis_successful.y.n. != "",
                            select = c(embryo_id, embryo_type, cytokinesis_successful.y.n.))

cytokinesis_phenotypes_grouped = cytokinesis_phenotypes %>%
  group_by(embryo_type, cytokinesis_successful.y.n.)%>%
  summarise(Count = n())

cytokinesis_sample_sizes = calculate_sample_sizes(cytokinesis_phenotypes, "cytokinesis_successful.y.n.")

cytokinesis_successful_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                            Hybrid"),"\n(N=", c(cytokinesis_sample_sizes$jk574,
                                                cytokinesis_sample_sizes$lkc28,
                                          cytokinesis_sample_sizes$hybrid),")",sep="")

embryo_descriptors = list("embryo_type","cytokinesis_successful.y.n.")

create_factors(cytokinesis_phenotypes_grouped, embryo_descriptors)

cytokinesis_phenotypes_visualization = ggplot(cytokinesis_phenotypes_grouped,
                                          mapping = aes(x = embryo_type, y = Count,
                                                        fill =cytokinesis_successful.y.n.))

 cytokinesis_phenotypes_visualization + geom_bar(position= "fill", stat= "identity") + xlab("Cross Type") + ylab("Fraction") + theme_cowplot(16) + scale_x_discrete(labels = cytokinesis_successful_Nlab)

ggsave(paste0(early_embryo_measurements_directory,"EHEP-DIC-GLOS-one_cell_cytokinesis_phenotypes.pdf"), width = 20, height = 20, units = "cm",
       device = "pdf", dpi = 300)



## ---------------------------------------------------------------------------------
cytoplasmic_clearings_phenotypes = one_cell_phenotypes[,c(1,2,12)]

cytoplasmic_clearings_phenotypes = subset(cytoplasmic_clearings_phenotypes, extra_cytoplasmic_clearings.y.n. != "",
                            select = c(embryo_id, embryo_type, extra_cytoplasmic_clearings.y.n.))

cytoplasmic_clearings_phenotypes_grouped = cytoplasmic_clearings_phenotypes %>%
  group_by(embryo_type, extra_cytoplasmic_clearings.y.n.)%>%
  summarise(Count = n())

embryo_descriptors = list("embryo_type","extra_cytoplasmic_clearings.y.n.")

create_factors(cytoplasmic_clearings_phenotypes_grouped, embryo_descriptors)

cytoplasmic_clearings_sample_sizes = calculate_sample_sizes(cytoplasmic_clearings_phenotypes, "extra_cytoplasmic_clearings.y.n.")

cytoplasmic_clearings_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                            Hybrid"),"\n(N=", c(cytoplasmic_clearings_sample_sizes$jk574,
                                                cytoplasmic_clearings_sample_sizes$lkc28,
                                          cytoplasmic_clearings_sample_sizes$hybrid),")",sep="")

cytoplasmic_clearings_phenotypes_visualization = ggplot(cytoplasmic_clearings_phenotypes_grouped,
                                          mapping = aes(x = embryo_type,
                                                        y = Count, fill = extra_cytoplasmic_clearings.y.n.))

cytoplasmic_clearings_phenotypes_visualization + geom_bar(position= "fill", stat= "identity") + xlab("Cross Type") + ylab("Fraction") + theme_cowplot(16) + scale_x_discrete(labels = cytoplasmic_clearings_Nlab)

ggsave(paste0(early_embryo_measurements_directory,"EHEP-DIC-GLOS-one_cell_cytoplasmic_clearings_phenotypes.pdf"), width = 20, height = 20, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------

source(paste0(function_directory,"select_variables_reformat_long.R"))
source(paste0(function_directory,"create_factors_sample_sizes.R"))
source(paste0(function_directory,"two_cell_size_comparisons.R"))
source(paste0(function_directory,"four_cell_phenotype_frequencies.R"))
source(paste0(function_directory,"select_variables_reformat_long.R"))
source(paste0(function_directory,"calculate_spindle_angle_ranges.R"))

spindle_angle_ranges = read.csv(paste0(early_embryo_measurements_directory,spindle_angle_range_measurements_filename))


## ---------------------------------------------------------------------------------
spindle_angle_ranges = drop_na(spindle_angle_ranges)

## convert dataframe to long 
spindle_angle_long = reshape(spindle_angle_ranges, 
        direction = "long",
        varying = c("spindle_angle","spindle_angle_range1","spindle_angle_range2","spindle_angle_range3", "spindle_angle_range4", "spindle_angle_range5"),
        v.names = "Angle_measurement",
        timevar = "spindle_measurement_time-point",
        times = c("spindle_angle","spindle_angle_range1","spindle_angle_range2","spindle_angle_range3", "spindle_angle_range4", "spindle_angle_range5"))

calculate_spindle_range = function(spindle_angle_measurements){
  spindle_angle_ranges = list()
  for (i in 1:nrow(spindle_angle_measurements)){
    print(i)
    angles = spindle_angle_measurements[i, 3:8]
    spindle_range = max(angles) - min(angles)
    spindle_angle_ranges = append(spindle_angle_ranges, spindle_range)
  }
  
  return(spindle_angle_ranges)
}
spindle_range_calculations = calculate_spindle_range(spindle_angle_ranges)
spindle_angle_ranges$spindle_ranges = spindle_range_calculations
spindle_angle_ranges$spindle_ranges = as.numeric(spindle_angle_ranges$spindle_ranges)

summary = spindle_angle_ranges %>% 
  group_by(EmbryoType) %>%
  summarise(mean = mean(spindle_ranges), se = sd(spindle_ranges)/sqrt(n()), n = n(), sd = sd(spindle_ranges))
print(summary)

spindle_angle_summary = spindle_angle_long %>%
  group_by(EmbryoType) %>%
  summarise(mean = mean(Angle_measurement), n= n(), sd = sd(Angle_measurement))
print(spindle_angle_summary)

#Calculate 2SD interval for wild-type embryos
wt_spindle_angles = subset(spindle_angle_long, EmbryoType != "hybrid")

wt_spindle_angles_avg = mean(wt_spindle_angles$Angle_measurement)
wt_spindle_angles_sd = sd(wt_spindle_angles$Angle_measurement)

wt_spindle_angles_conf_int_upper = wt_spindle_angles_avg + 2*wt_spindle_angles_sd
wt_spindle_angles_conf_int_lower = wt_spindle_angles_avg - 2*wt_spindle_angles_sd

#calculate wt spindle angles that fall in wt 95% data interval
wt_spindle_angles_in_wt_range = subset(wt_spindle_angles,
                                           wt_spindle_angles$Angle_measurement >
                                             wt_spindle_angles_conf_int_lower &
                                             wt_spindle_angles$Angle_measurement <
                                             wt_spindle_angles_conf_int_upper)

# calculate hybrid spindle angles that fall into wt 95% interval

hybrid_spindle_angles = subset(spindle_angle_long, spindle_angle_long$EmbryoType == "hybrid")
hybrid_spindle_angles_in_wt_range = subset(hybrid_spindle_angles,
                                           hybrid_spindle_angles$Angle_measurement >
                                             wt_spindle_angles_conf_int_lower &
                                             hybrid_spindle_angles$Angle_measurement <
                                             wt_spindle_angles_conf_int_upper)
hybrid_spindle_angles_n = nrow(hybrid_spindle_angles)
hybrid_spindle_angles_in_wt_range_n = nrow(hybrid_spindle_angles_in_wt_range)
percent_hybrids_in_95int_wt = hybrid_spindle_angles_in_wt_range_n/hybrid_spindle_angles_n

spindle_angle_range_Nlab = paste(c("hybrid",
                                "c. brenneri",
                                "c. elegans"),"\n(N=", c(summary$n[summary$EmbryoType == "hybrid"],
                                                summary$n[summary$EmbryoType == "C. brenneri"],
                                          summary$n[summary$EmbryoType == "C. elegans"]),")",sep="")


## ---------------------------------------------------------------------------------
hybrid_spindle_angle_range = range(hybrid_spindle_angles$Angle_measurement)


## ---------------------------------------------------------------------------------
#1. separate data into hybrid, cbn, and cel dataframes
hybrid_spindle_ranges = subset(spindle_angle_ranges, EmbryoType == "hybrid")
cbn_spindle_ranges = subset(spindle_angle_ranges, EmbryoType == "C. brenneri")
cel_spindle_ranges = subset(spindle_angle_ranges, EmbryoType == "C. elegans")
#2. count unique instances of names in each subset
number_of_hyb_embryos = length(unique(hybrid_spindle_ranges$EmbryoID))
number_of_cbn_embryos = length(unique(cbn_spindle_ranges$EmbryoID))
number_of_cel_embryos = length(unique(cel_spindle_ranges$EmbryoID))


## ---------------------------------------------------------------------------------
spindle_angle_ranges$EmbryoType = factor(spindle_angle_ranges$EmbryoType,
                                                             levels= unique(spindle_angle_ranges$EmbryoType))

spindle_angle_ranges = arrange(spindle_angle_ranges, EmbryoType)
spindle_angle_long_arranged = arrange(spindle_angle_long, EmbryoType)
spindle_angle_long_arranged_names = unique(spindle_angle_long_arranged$EmbryoID)


## ---------------------------------------------------------------------------------
spindle_angle_range_vis = ggplot(spindle_angle_long_arranged, aes(EmbryoID, Angle_measurement, color = EmbryoType))

spindle_angle_range_vis + geom_point(alpha = 0.6, size = 3) + theme_cowplot(14) +   scale_color_manual(values = c("hybrid" = "purple2", "C. brenneri"="black","C. elegans"="grey59")) + geom_hline(yintercept = 0, col = "blue", linetype='dotted') + scale_x_discrete(limits = spindle_angle_long_arranged_names)

## ---------------------------------------------------------------------------------
ggsave(paste0(early_embryo_measurements_directory,"EHEP-DIC-GLOS-individual_spindle_angles.pdf"), width = 20, height = 20, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
spindle_range_vis = ggplot(summary, aes(EmbryoType, mean))

spindle_range_vis + geom_point() + geom_errorbar(aes(ymin=mean-se,
                              ymax =mean + se), width = 0.2,
                                    position = position_dodge(0.05)) + theme_cowplot(16) + scale_x_discrete(labels= spindle_angle_range_Nlab)

## ---------------------------------------------------------------------------------
ggsave(paste0(early_embryo_measurements_directory,"EHEP-DIC-GLOS-spindle_angle_ranges.pdf"), width = 20, height = 20, units = "cm",
       device = "pdf", dpi = 300)

