## Visualize quantified phenotypes from IF experiments

library(ggplot2); library(dplyr); library(reshape); library(Hmisc); library(tidyverse);
library(stats); library(EnvStats); library(cowplot)

get_file_info = function(trailingOnly = TRUE, asValues = TRUE) {
  data_file_info = commandArgs(TRUE)
  if (length(data_file_info) < 3){
    stop("Provide the file path to functions,
         data file dir, immunofluorescence measurements,
         eg. dev/ IF_measurements_data_dir/ IF_quantifications.ext")
  } else {
    return(data_file_info)
  }
}

immunofluorescence_measurement_file_info = get_file_info()

function_directory = immunofluorescence_measurement_file_info[1]
early_embryo_measurements_directory = immunofluorescence_measurement_file_info[2]
immunofluorescence_measurements_filename = immunofluorescence_measurement_file_info[3]

source(paste0(function_directory,"select_variables_reformat_long.R"))
source(paste0(function_directory,"create_factors_sample_sizes.R"))
source(paste0(function_directory,"two_cell_size_comparisons.R"))
source(paste0(function_directory,"four_cell_phenotype_frequencies.R"))
source(paste0(function_directory,"select_variables_reformat_long.R"))

## import data
immunofluorescence_measurements = read.csv(paste0(early_embryo_measurements_directory,
                                                 immunofluorescence_measurements_filename), header = TRUE)

meiotic_failure_and_capture_measurements = immunofluorescence_measurements[,c(1,5,6)]
meiotic_failure_and_capture_measurements = drop_na(meiotic_failure_and_capture_measurements)
meiotic_failure_measurements = immunofluorescence_measurements[,c(1,6)]
meiotic_failure_measurements = drop_na(meiotic_failure_measurements)

spindle_arm_quantifications = immunofluorescence_measurements[, c(1,2,7)]
spindle_arm_quantifications = drop_na(spindle_arm_quantifications)

spindle_arm_and_meiotic_failure = immunofluorescence_measurements[, c(1,2,6,7)]
spindle_arm_and_meiotic_failure = drop_na(spindle_arm_and_meiotic_failure)

## group variables for meiotic capture type

grouped_meiotic_failure_and_capture = meiotic_failure_and_capture_measurements %>%
  group_by(meiotic_DNA_capture, meiotic_failure)%>%
  summarise(all_embs_count=n())

grouped_meiotic_failure_measurements = meiotic_failure_measurements %>%
  group_by(meiotic_failure) %>%
  summarise(all_embs_count=n())

grouped_spindle_arm_quantifications = spindle_arm_quantifications %>%
  group_by(spindle_arms, embryo_type) %>%
  summarise(embryo_number = n())

grouped_spindle_arm_and_meiotic_failure = spindle_arm_and_meiotic_failure %>%
  group_by(meiotic_failure, spindle_arms, embryo_type) %>%
  summarise(embryo_number = n())

## create factors for visualization
data_descriptors = list("meiotic_failure", "meiotic_DNA_capture")
create_factors(grouped_meiotic_failure_and_capture, data_descriptors = data_descriptors)

meiotic_failure_data_descriptors = list("meiotic_failure")
create_factors(grouped_meiotic_failure_measurements, data_descriptors = meiotic_failure_data_descriptors)

meiotic_failure_and_dna_capture_n = c(grouped_meiotic_failure_and_capture$all_embs_count[1] +
                                        grouped_meiotic_failure_and_capture$all_embs_count[4],
                                      grouped_meiotic_failure_and_capture$all_embs_count[2],
                                      grouped_meiotic_failure_and_capture$all_embs_count[3])

dna_capture_n = sum(subset(grouped_meiotic_failure_and_capture,
                           grouped_meiotic_failure_and_capture$meiotic_DNA_capture == "y")$all_embs_count)
dna_not_captured_n = sum(subset(grouped_meiotic_failure_and_capture,
                                grouped_meiotic_failure_and_capture$meiotic_DNA_capture == "n")$all_embs_count)

spindle_data_descriptors = list("spindle_arms")
create_factors(grouped_spindle_arm_quantifications, data_descriptors = spindle_data_descriptors)

spindle_arms_n = subset(grouped_spindle_arm_quantifications,
                        grouped_spindle_arm_quantifications$spindle_arms == "y")$embryo_number
no_spindle_arms_n = subset(grouped_spindle_arm_quantifications,
                           grouped_spindle_arm_quantifications$spindle_arms == "n")$embryo_number


## define sample sizes
meiotic_failure_and_capture_n = paste(c("meiosis I failure",
                             "meiosis II failure",
                             "no failure"),"\n(N=", c(meiotic_failure_and_dna_capture_n[3],
                                                       meiotic_failure_and_dna_capture_n[1],
                                                       meiotic_failure_and_dna_capture_n[2]),")",sep="")

meiotic_failure_n = paste(c("meiosis I failure",
                            "meiosis II failure",
                            "no failure"),"\n(N=", c(grouped_meiotic_failure_measurements$all_embs_count[3],
                                                     grouped_meiotic_failure_measurements$all_embs_count[1],
                                                     grouped_meiotic_failure_measurements$all_embs_count[2]),")",sep="")

dna_capture_no_capture_n = paste(c("no dna capture", 
                                    "dna capture"), "\n(N=", 
                                  c(dna_not_captured_n,
                                    dna_capture_n),")", sep="")

spindle_arm_n_label = paste(c("no spindle arms", 
                              "spindle arms"), "\n(N=", 
                            c(no_spindle_arms_n,
                              spindle_arms_n),")", sep="")
## visualize meiotic spindle failure
meiotic_failure_and_dna_capture_vis = ggplot(grouped_meiotic_failure_and_capture, mapping = aes(meiotic_failure,
                                                                                                all_embs_count,
                                                                                                fill = meiotic_DNA_capture))

meiotic_failure_and_dna_capture_vis + geom_bar(position= "stack", stat= "identity") + theme_cowplot(14) +
  scale_x_discrete(labels= meiotic_failure_and_capture_n)

ggsave(paste0(early_embryo_measurements_directory,"meiotic_dna_failure_capture.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)

meiotic_failure_vis = ggplot(grouped_meiotic_failure_measurements, mapping = aes(meiotic_failure,
                                                                                                all_embs_count))

meiotic_failure_vis + geom_bar(position= "dodge", stat= "identity") + theme_cowplot(14) +
  scale_x_discrete(labels= meiotic_failure_n)

ggsave(paste0(early_embryo_measurements_directory,"meiotic_dna_failure.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


dna_capture_vis = ggplot(grouped_meiotic_failure_and_capture, mapping = aes(meiotic_DNA_capture, 
                                                                            all_embs_count, 
                                                                            fill = meiotic_failure))
dna_capture_vis + geom_bar(position = "fill", stat = "identity") + theme_cowplot(14) +
  scale_x_discrete(labels = dna_capture_no_capture_n)

ggsave(paste0(early_embryo_measurements_directory,"dna_capture_and_meiotic_failure.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)

spindle_arm_vis = ggplot(grouped_spindle_arm_quantifications, mapping = aes(embryo_type,
                                                                            embryo_number,
                                                                            fill = spindle_arms))
spindle_arm_vis + geom_bar(position = "stack", stat = "identity") + theme_cowplot(14) + 
  scale_x_discrete(labels = spindle_arm_n_label)

spindle_arm_vis + geom_bar(position = "fill", stat = "identity") + theme_cowplot(14) + 
  scale_x_discrete(labels = spindle_arm_n_label)

ggsave(paste0(early_embryo_measurements_directory,"spindle_arms_stacked.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)