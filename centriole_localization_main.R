## centriole localization main

## Import libraries

library(ggplot2); library(dplyr); library(reshape); library(Hmisc);
library(tidyverse); library(stats); library(EnvStats); library(cowplot); library(glue)

## Import global user functions
get_file_info = function(trailingOnly = TRUE, asValues = TRUE) {
  data_file_info = commandArgs(TRUE)
  if (length(data_file_info) < 3){
    stop("Provide the file path to functions,
         data file dir, and filename, eg. dev/ centriole_data_dir/ centriole_data.ext")
  } else {
    return(data_file_info)
  }
}

centriole_file_info = get_file_info()
setwd(centriole_file_info[1])
source("get_functions.R")
user_functions = list("select_variables_reformat_long.R","create_factors_sample_sizes.R",
                      "centriole_localization_data_handling.R")

get_functions(centriole_file_info[1], user_functions)

centriole_file_directory = centriole_file_info[2]
centriole_filename = centriole_file_info[3]
centriole_file_path = paste0(centriole_file_directory,centriole_filename)
## Import data
centriole_qualification = read.csv(centriole_file_path)
broad_position = centriole_qualification[, -c(4, 5, 6)]
broad_position = broad_position[!grepl("06723", broad_position$emb_id), ]
broad_position = broad_position[!grepl("12723", broad_position$emb_id), ]
broad_position = broad_position[!grepl("02823", broad_position$emb_id), ]
broad_position = broad_position[!grepl("071223", broad_position$emb_id), ]
broad_position = subset(broad_position, emb_type != "brenneri")
broad_position_descriptors = c("emb_type")
create_factors(broad_position, broad_position_descriptors)
broad_position = broad_position %>%
  group_by(emb_type, aberrant_positioning) %>%
  summarise(Count = n())

centrosome_errors = centriole_qualification[, -c(3,5,6)]
centrosome_errors = centrosome_errors[!grepl("06723", centrosome_errors$emb_id), ]
centrosome_errors = centrosome_errors[!grepl("12723", centrosome_errors$emb_id), ]
centrosome_errors = centrosome_errors[!grepl("02823", centrosome_errors$emb_id), ]
centrosome_errors = centrosome_errors[!grepl("071223", centrosome_errors$emb_id), ]
centrosome_errors = subset(centrosome_errors, emb_type != "brenneri")
centrosome_errors_descriptors = c("emb_type")
create_factors(centrosome_errors, centrosome_errors_descriptors)
centrosome_errors = centrosome_errors %>%
  group_by(emb_type, error1) %>%
  summarise(Count = n())

centrosome_detachment = centriole_qualification[ , c(1:2,7)]
centrosome_detachment = centrosome_detachment[!grepl("06723", centrosome_detachment$emb_id), ]
centrosome_detachment = centrosome_detachment[!grepl("12723", centrosome_detachment$emb_id), ]
centrosome_detachment = centrosome_detachment[!grepl("02823", centrosome_detachment$emb_id), ]
centrosome_detachment = centrosome_detachment[!grepl("071223", centrosome_detachment$emb_id), ]
centrosome_detachment = subset(centrosome_detachment, emb_type != "brenneri")
centrosome_detachment_descriptors = c("emb_type")
create_factors(centrosome_detachment, centrosome_detachment_descriptors)
centrosome_detachment_summary = centrosome_detachment %>%
  group_by(emb_type, detachment_before_NEBD) %>%
  summarise(Count = n())

centrosome_nuclei_position = centriole_qualification[ , c(1:2,8)]
centrosome_nuclei_position = centrosome_nuclei_position[!grepl("06723", centrosome_nuclei_position$emb_id), ]
centrosome_nuclei_position = centrosome_nuclei_position[!grepl("12723", centrosome_nuclei_position$emb_id), ]
centrosome_nuclei_position = centrosome_nuclei_position[!grepl("02823", centrosome_nuclei_position$emb_id), ]
centrosome_nuclei_position = centrosome_nuclei_position[!grepl("071223", centrosome_nuclei_position$emb_id), ]
centrosome_nuclei_position = subset(centrosome_nuclei_position, emb_type != "brenneri")
centrosome_nuclei_position_descriptors = c("emb_type")
create_factors(centrosome_nuclei_position, centrosome_nuclei_position_descriptors)
centrosome_nuclei_position_summary = centrosome_nuclei_position %>%
  group_by(emb_type, between_pronuclei_at_pn_meeting) %>%
  summarise(Count = n())

## segregate nuclear size data
pronuclear_size = centriole_qualification[ , c(1:2,7,10:15)]
pronuclear_size = pronuclear_size[!grepl("06723", pronuclear_size$emb_id), ]
pronuclear_size = pronuclear_size[!grepl("12723", pronuclear_size$emb_id), ]
pronuclear_size = pronuclear_size[!grepl("02823", pronuclear_size$emb_id), ]
pronuclear_size = pronuclear_size[!grepl("071223", pronuclear_size$emb_id), ]

pronuclear_size_pn_meeting = pronuclear_size[, c(1:3, 4:7)]
pronuclear_size_pn_meeting = drop_na(pronuclear_size_pn_meeting)
pronuclear_size_before_pn_meeting = pronuclear_size[, c(1:3, 8:9)]
pronuclear_size_before_pn_meeting = drop_na(pronuclear_size_before_pn_meeting)
## manipulate data
female_pn_area = pi*(pronuclear_size_pn_meeting$female_pn_height/2)*(pronuclear_size_pn_meeting$female_pn_length/2)
male_pn_area = pi*(pronuclear_size_pn_meeting$male_pn_height/2)*(pronuclear_size_pn_meeting$male_pn_length/2)

male_pn_before_pn_meeting_area = pi*(pronuclear_size_before_pn_meeting$male_pn_height_before_pn_meeting/2)*(pronuclear_size_before_pn_meeting$male_pn_length_before_pn_meeting/2)

pronuclear_size_pn_meeting$female_pn_area = female_pn_area
pronuclear_size_pn_meeting$male_pn_area = male_pn_area
pronuclear_size_pn_meeting$male_to_female_pn_area = pronuclear_size_pn_meeting$male_pn_area/pronuclear_size_pn_meeting$female_pn_area

pronuclear_size_before_pn_meeting$early_pn_area = male_pn_before_pn_meeting_area


## convert both dataframes to long
pn_area_pn_meeting_long = reshape(pronuclear_size_pn_meeting, 
                             direction = "long",
                             varying = c("female_pn_height","female_pn_length","male_pn_height","male_pn_length", "female_pn_area", "male_pn_area", "male_to_female_pn_area"),
                             v.names = "pn_size_measurement",
                             timevar = "measurement_type",
                             times = c("female_pn_height","female_pn_length","male_pn_height","male_pn_length", "female_pn_area", "male_pn_area", "male_to_female_pn_area"))

pn_area_before_pn_meeting_long = reshape(pronuclear_size_before_pn_meeting, 
                                  direction = "long",
                                  varying = c("male_pn_height_before_pn_meeting","male_pn_length_before_pn_meeting","early_pn_area"),
                                  v.names = "pn_size_measurement",
                                  timevar = "measurement_type",
                                  times = c("male_pn_height_before_pn_meeting","male_pn_length_before_pn_meeting","early_pn_area"))

## calculate median early pn area
cel_early_pn_area_median = median(subset(pronuclear_size_before_pn_meeting, emb_type == "elegans")$early_pn_area)
hybrid_early_pn_area_median = median(subset(pronuclear_size_before_pn_meeting, emb_type == "hybrid")$early_pn_area)
cbn_early_pn_area_median = median(subset(pronuclear_size_before_pn_meeting, emb_type == "brenneri")$early_pn_area)

## calculate IQR for early pn area
cel_early_pn_area_IQR = IQR(subset(pronuclear_size_before_pn_meeting, emb_type == "elegans")$early_pn_area)
cbn_early_pn_area_IQR = IQR(subset(pronuclear_size_before_pn_meeting, emb_type == "brenneri")$early_pn_area)
hybrid_early_pn_area_IQR = IQR(subset(pronuclear_size_before_pn_meeting, emb_type == "hybrid")$early_pn_area)

## calculate median and IQR for ratio of sperm pn area to oocyte pn area
cel_pn_area_meeting_median = median(subset(pronuclear_size_pn_meeting, emb_type == "elegans")$male_to_female_pn_area)
hybrid_pn_area_meeting_median = median(subset(pronuclear_size_pn_meeting, emb_type == "hybrid")$male_to_female_pn_area)
cbn_pn_area_meeting_median = median(subset(pronuclear_size_pn_meeting, emb_type == "brenneri")$male_to_female_pn_area)

cel_pn_area_meeting_IQR = IQR(subset(pronuclear_size_pn_meeting, emb_type == "elegans")$male_to_female_pn_area)
hybrid_pn_area_meeting_IQR = IQR(subset(pronuclear_size_pn_meeting, emb_type == "hybrid")$male_to_female_pn_area)
cbn_pn_area_meeting_IQR = IQR(subset(pronuclear_size_pn_meeting, emb_type == "brenneri")$male_to_female_pn_area)

## calculate median and IQR for female and male pn area at pn meeting
cel_pn_area_meeting_female_median = median(subset(pronuclear_size_pn_meeting, emb_type == "elegans")$female_pn_area)
hybrid_pn_area_meeting_female_median = median(subset(pronuclear_size_pn_meeting, emb_type == "hybrid")$female_pn_area)
cbn_pn_area_meeting_female_median = median(subset(pronuclear_size_pn_meeting, emb_type == "brenneri")$female_pn_area)

cel_pn_area_meeting_female_IQR = IQR(subset(pronuclear_size_pn_meeting, emb_type == "elegans")$female_pn_area)
hybrid_pn_area_meeting_female_IQR = IQR(subset(pronuclear_size_pn_meeting, emb_type == "hybrid")$female_pn_area)
cbn_pn_area_meeting_female_IQR = IQR(subset(pronuclear_size_pn_meeting, emb_type == "brenneri")$female_pn_area)

cel_pn_area_meeting_male_median = median(subset(pronuclear_size_pn_meeting, emb_type == "elegans")$male_pn_area)
hybrid_pn_area_meeting_male_median = median(subset(pronuclear_size_pn_meeting, emb_type == "hybrid")$male_pn_area)
cbn_pn_area_meeting_male_median = median(subset(pronuclear_size_pn_meeting, emb_type == "brenneri")$male_pn_area)

cel_pn_area_meeting_male_IQR = IQR(subset(pronuclear_size_pn_meeting, emb_type == "elegans")$male_pn_area)
hybrid_pn_area_meeting_male_IQR = IQR(subset(pronuclear_size_pn_meeting, emb_type == "hybrid")$male_pn_area)
cbn_pn_area_meeting_male_IQR = IQR(subset(pronuclear_size_pn_meeting, emb_type == "brenneri")$male_pn_area)




wt_male_pronuclear_size_pn_meeting = subset(pn_area_pn_meeting_long, emb_type != "hybrid" & measurement_type == "male_pn_area")
avg_wt_male_pronuclear_size_pn_meeting = mean(wt_male_pronuclear_size_pn_meeting$pn_size_measurement)
sd_wt_male_pronuclear_size_pn_meeting = sd(wt_male_pronuclear_size_pn_meeting$pn_size_measurement)
upper_conf_int = avg_wt_male_pronuclear_size_pn_meeting + 2*sd_wt_male_pronuclear_size_pn_meeting
lower_conf_int = avg_wt_male_pronuclear_size_pn_meeting - 2*sd_wt_male_pronuclear_size_pn_meeting

cbn_male_pronuclear_size_pn_meeting = subset(pn_area_pn_meeting_long, emb_type == "brenneri" & measurement_type == "male_pn_area")
avg_cbn_male_pronuclear_size_pn_meeting = mean(cbn_male_pronuclear_size_pn_meeting$pn_size_measurement)
sd_cbn_male_pronuclear_size_pn_meeting = sd(cbn_male_pronuclear_size_pn_meeting$pn_size_measurement)
cbn_upper_conf_int = avg_cbn_male_pronuclear_size_pn_meeting + 2*sd_cbn_male_pronuclear_size_pn_meeting
cbn_lower_conf_int = avg_cbn_male_pronuclear_size_pn_meeting - 2*sd_cbn_male_pronuclear_size_pn_meeting


hybrid_male_pronuclear_size_pn_meeting = subset(pn_area_pn_meeting_long, emb_type == "hybrid" & measurement_type == "male_pn_area")
hybrid_male_pn_in_wt_95int = subset(hybrid_male_pronuclear_size_pn_meeting, pn_size_measurement > lower_conf_int & pn_size_measurement < upper_conf_int)
no_of_hybrid_male_pn_in_wt_95int = nrow(hybrid_male_pn_in_wt_95int)

hybrid_male_pn_in_cbn_95int = subset(hybrid_male_pronuclear_size_pn_meeting, pn_size_measurement > cbn_lower_conf_int
                                     & pn_size_measurement < cbn_upper_conf_int)
# create factors for graphing detachment on pn size
pn_area_pn_meeting_descriptors = list("emb_type", "detachment_before_NEBD", "measurement_type")
create_factors(pn_area_pn_meeting_long, pn_area_pn_meeting_descriptors)
pronuclear_meeting_elegans_N = sum(pronuclear_size_pn_meeting$emb_type == "elegans")
pronuclear_meeting_hybrid_N = sum(pronuclear_size_pn_meeting$emb_type == "hybrid")
pronuclear_meeting_brenneri_N = sum(pronuclear_size_pn_meeting$emb_type == "brenneri")
pn_meeting_area_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                                Hybrid"),
                              "\n(N=",c(pronuclear_meeting_elegans_N,
                                        pronuclear_meeting_hybrid_N,
                                        pronuclear_meeting_brenneri_N)) 

pn_area_before_pn_meeting_descriptors = list("emb_type", "detachment_before_NEBD", "measurement_type")
create_factors(pn_area_before_pn_meeting_long, pn_area_before_pn_meeting_descriptors)

pronuclear_before_pn_meeting_elegans_N = sum(pronuclear_size_before_pn_meeting$emb_type == "elegans")
pronuclear_before_pn_meeting_hybrid_N = sum(pronuclear_size_before_pn_meeting$emb_type == "hybrid")
pronuclear_before_pn_meeting_brenneri_N = sum(pronuclear_size_before_pn_meeting$emb_type == "brenneri")

pn_before_meeting_area_Nlab = paste(c("C. elegans Wild-type",
                               "C. brenneri(f) x C. elegans(m) 
                                Hybrid"),
                             "\n(N=",c(pronuclear_before_pn_meeting_elegans_N,
                                       pronuclear_before_pn_meeting_hybrid_N,
                                       pronuclear_before_pn_meeting_brenneri_N)) 

## centriole segregation
centrosome_segregation = centriole_qualification[ , c(1:2,16)]
centrosome_segregation = centrosome_segregation[!grepl("06723", centrosome_segregation$emb_id), ]
centrosome_segregation = centrosome_segregation[!grepl("12723", centrosome_segregation$emb_id), ]
centrosome_segregation = centrosome_segregation[!grepl("02823", centrosome_segregation$emb_id), ]
centrosome_segregation = centrosome_segregation[!grepl("071223", centrosome_segregation$emb_id), ]
centrosome_segregation = subset(centrosome_segregation, emb_type != "brenneri")
centrosome_segregation_descriptors = c("emb_type")
create_factors(centrosome_segregation, centrosome_segregation_descriptors)
centrosome_segregation_summary = centrosome_segregation %>%
  group_by(emb_type, centrioles_distribute_two_cells) %>%
  summarise(Count = n())
# visualize data
broad_position_vis = ggplot(broad_position, mapping = aes(emb_type, Count, fill = aberrant_positioning))
broad_position_vis + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14)
ggsave(paste0(centriole_file_info[2],"broad_centrosome_positioning.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)

centrosome_error_vis = ggplot(centrosome_errors, mapping = aes(emb_type, Count, fill = error1))
centrosome_error_vis + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14)
ggsave(paste0(centriole_file_info[2],"centrosome_error.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)

centrosome_detachment_vis = ggplot(centrosome_detachment_summary, mapping = aes(emb_type, Count, fill = detachment_before_NEBD))
centrosome_detachment_vis + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14)
ggsave(paste0(centriole_file_info[2],"centrosome_detachment.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)

nuclei_positioning_vis = ggplot(centrosome_nuclei_position_summary, mapping = aes(emb_type, Count, fill = between_pronuclei_at_pn_meeting))
nuclei_positioning_vis + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14)
ggsave(paste0(centriole_file_info[2],"centrosome_positioning_at_pn_meeting.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)

pn_area_pn_meeting_vis = ggplot(subset(pn_area_pn_meeting_long, measurement_type == "male_pn_area"|measurement_type == "female_pn_area"),
                                aes(emb_type, pn_size_measurement))
pn_area_pn_meeting_vis + geom_boxplot() +
  geom_point(position = position_jitter(w=0.3, h=0), size = 5, alpha = 0.7, aes(color = detachment_before_NEBD)) + facet_wrap(~measurement_type) +
  theme_cowplot(14) + scale_y_continuous(limits = c(0, 85)) + labs(y = "pn area (uM2)", x = "Cross") +
  scale_fill_manual(values = c("black", "grey", "purple")) + scale_x_discrete(labels= pn_meeting_area_Nlab)

ggsave(paste0(centriole_file_info[2],"pronuclear_size_at_pn_meeting_detachment_overlay.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)

pn_area_pn_meeting_ratio_vis = ggplot(subset(pn_area_pn_meeting_long, measurement_type == "male_to_female_pn_area"),
                                      aes(emb_type, pn_size_measurement))
pn_area_pn_meeting_ratio_vis + geom_boxplot() +
  geom_point(position = position_jitter(w=0.2, h=0), size = 5, alpha = 0.6) +
  theme_cowplot(14) + scale_y_continuous(limits = c(0, 2)) + labs(y = "pn area ratio", x = "Cross") +
  scale_color_grey() + scale_x_discrete(labels= pn_meeting_area_Nlab)

ggsave(paste0(centriole_file_info[2],"pronuclear_size_at_pn_meeting_ratio.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)

pn_area_before_pn_meeting_vis = ggplot(subset(pn_area_before_pn_meeting_long, measurement_type == "early_pn_area"),
                                       aes(emb_type, pn_size_measurement))
pn_area_before_pn_meeting_vis + geom_boxplot() +
  geom_point(position = position_jitter(w=0.2, h=0), size = 5, alpha = 0.7, aes(color = detachment_before_NEBD)) +
  theme_cowplot(14) + scale_y_continuous(limits = c(0, 70)) + labs(y = "pn area (uM2)", x = "Cross") +
  scale_color_grey() + scale_x_discrete(labels= pn_before_meeting_area_Nlab)

ggsave(paste0(centriole_file_info[2],"pronuclear_size_before_pn_meeting_detachment_overlays.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)

centrosome_segregation_vis = ggplot(centrosome_segregation_summary, mapping = aes(emb_type, Count, fill = centrioles_distribute_two_cells))
centrosome_segregation_vis + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_fill_grey()
ggsave(paste0(centriole_file_info[2],"centrosome_positioning_two_cell.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)
