## ----setup, include=FALSE-----------------------
knitr::opts_chunk$set(echo = TRUE)


## -----------------------------------------------
library(ggplot2); library(dplyr); library(reshape); library(Hmisc);
library(tidyverse); library(stats); library(EnvStats); library(cowplot)

get_file_info = function(trailingOnly = TRUE, asValues = TRUE) {
  data_file_info = commandArgs(TRUE)
  if (length(data_file_info) < 4){
    stop("Provide the file path to functions,
         data file dir, and filename, eg. dev/ blinded_embryos_scores_data_dir/ blinded_embryo_scores_data.ext blinded_embryo_scores_for_heat_map_from_python.ext")
  } else {
    return(data_file_info)
  }
}

blinded_embryo_scores_file_info = get_file_info()
function_directory = blinded_embryo_scores_file_info[1]
blinded_embryo_scores_file_directory = blinded_embryo_scores_file_info[2]
blinded_embryo_scores_filename = blinded_embryo_scores_file_info[3]
blinded_embryo_scores_for_heat_map = blinded_embryo_scores_file_info[4]

## Import global user functions
setwd(function_directory)
source("get_functions.R")
user_functions = list("select_variables_reformat_long.R","create_factors_sample_sizes.R",
                      "centriole_localization_data_handling.R")

get_functions(function_directory, user_functions)


## -----------------------------------------------
embryo_scores = read.csv(paste0(blinded_embryo_scores_file_directory,blinded_embryo_scores_filename))


## -----------------------------------------------
embryo_scores[embryo_scores == "cant assess"] = NA
embryo_scores = filter(embryo_scores, embryo_type != "cni_hybrid" & embryo_type != "ju1325")


## -----------------------------------------------
polar_body_size = embryo_scores[, c(2, 3, 4)]
polar_body_size = drop_na(polar_body_size)
polar_body_size_grouped =  polar_body_size %>%
  group_by(embryo_type, irregular_polarbody) %>% summarise(Count = n())


## -----------------------------------------------
polar_body_sample_sizes = calculate_sample_sizes(polar_body_size, "irregular_polarbody")
sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(polar_body_sample_sizes$c48_hybrid,
                           polar_body_sample_sizes$cel_hybrid,
                           polar_body_sample_sizes$cre_hybrid,
                           polar_body_sample_sizes$csp48,
                          polar_body_sample_sizes$em464,
                          polar_body_sample_sizes$jk574, polar_body_sample_sizes$lkc28),")",sep="")

polar_body_descriptors = list("embryo_type")
create_factors(polar_body_size_grouped, polar_body_descriptors)


## -----------------------------------------------
polar_body_visualization = ggplot(polar_body_size_grouped,
                                mapping = aes(x = embryo_type, y = Count, fill = irregular_polarbody))

polar_body_visualization + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)


## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_polar_body.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
pronuclear_size = embryo_scores[, c(2,3,5)]
pronuclear_size = drop_na(pronuclear_size)
pronuclear_size_grouped =  pronuclear_size %>%
  group_by(embryo_type, uneven_sized_nuclei_pronuclear.meeting) %>% summarise(Count = n())

## -----------------------------------------------
pn_size_sample_sizes = calculate_sample_sizes(pronuclear_size, "uneven_sized_nuclei_pronuclear.meeting")
sample_size_Nlab = sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(pn_size_sample_sizes$c48_hybrid,
                           pn_size_sample_sizes$cel_hybrid,
                           pn_size_sample_sizes$cre_hybrid,
                           pn_size_sample_sizes$csp48,
                          pn_size_sample_sizes$em464,
                          pn_size_sample_sizes$jk574, pn_size_sample_sizes$lkc28),")",sep="")

pn_size_descriptors = list("embryo_type")
create_factors(pronuclear_size_grouped, polar_body_descriptors)


## -----------------------------------------------
pn_size_visualization = ggplot(pronuclear_size_grouped,
                                mapping = aes(x = embryo_type, y = Count, fill = uneven_sized_nuclei_pronuclear.meeting))

pn_size_visualization + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)



## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_pn_sizes.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
pronuclear_meeting = embryo_scores[, c(2,3,6)]
pronuclear_meeting = drop_na(pronuclear_meeting)
pronuclear_meeting_grouped =  pronuclear_meeting %>%
  group_by(embryo_type, pronuclei_fail_to_meet) %>% summarise(Count = n())


## -----------------------------------------------
pn_meeting_sample_sizes = calculate_sample_sizes(pronuclear_meeting, "pronuclei_fail_to_meet")
sample_size_Nlab = sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(pn_meeting_sample_sizes$c48_hybrid,
                           pn_meeting_sample_sizes$cel_hybrid,
                           pn_meeting_sample_sizes$cre_hybrid,
                           pn_meeting_sample_sizes$csp48,
                          pn_meeting_sample_sizes$em464,
                          pn_meeting_sample_sizes$jk574, pn_meeting_sample_sizes$lkc28),")",sep="")

pn_meeting_descriptors = list("embryo_type")
create_factors(pronuclear_meeting_grouped, pn_meeting_descriptors)


## -----------------------------------------------
pn_meeting_visualization = ggplot(pronuclear_meeting_grouped,
                                mapping = aes(x = embryo_type, y = Count, fill = pronuclei_fail_to_meet))

pn_meeting_visualization + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)



## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_pn_meeting.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
pronuclear_shape = embryo_scores[, c(2,3,7)]
pronuclear_shape = drop_na(pronuclear_shape)
pronuclear_shape_grouped =  pronuclear_shape %>%
  group_by(embryo_type, misshappen_nuclei) %>% summarise(Count = n())


## -----------------------------------------------
pn_shape_sample_sizes = calculate_sample_sizes(pronuclear_shape, "misshappen_nuclei")
sample_size_Nlab = sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(pn_shape_sample_sizes$c48_hybrid,
                           pn_shape_sample_sizes$cel_hybrid,
                           pn_shape_sample_sizes$cre_hybrid,
                           pn_shape_sample_sizes$csp48,
                          pn_shape_sample_sizes$em464,
                          pn_shape_sample_sizes$jk574, pn_shape_sample_sizes$lkc28),")",sep="")

pn_shape_descriptors = list("embryo_type")
create_factors(pronuclear_shape_grouped, pn_shape_descriptors)


## -----------------------------------------------
pn_shape_visualization = ggplot(pronuclear_shape_grouped,
                                mapping = aes(x = embryo_type, y = Count, fill = misshappen_nuclei))

pn_shape_visualization + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)



## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_pn_shape.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
pronuclear_movement = embryo_scores[, c(2,3,8)]
pronuclear_movement = drop_na(pronuclear_movement)
pronuclear_movement_grouped =  pronuclear_movement %>%
  group_by(embryo_type, abnormal_pronuclei_movement) %>% summarise(Count = n())


## -----------------------------------------------
pn_movement_sample_sizes = calculate_sample_sizes(pronuclear_movement, "abnormal_pronuclei_movement")
sample_size_Nlab = sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(polar_body_sample_sizes$c48_hybrid,
                           polar_body_sample_sizes$cel_hybrid,
                           polar_body_sample_sizes$cre_hybrid,
                           polar_body_sample_sizes$csp48,
                          polar_body_sample_sizes$em464,
                          polar_body_sample_sizes$jk574, polar_body_sample_sizes$lkc28),")",sep="")

pn_movement_descriptors = list("embryo_type")
create_factors(pronuclear_movement_grouped, pn_movement_descriptors)


## -----------------------------------------------
pn_movement_visualization = ggplot(pronuclear_movement_grouped,
                                mapping = aes(x = embryo_type, y = Count, fill = abnormal_pronuclei_movement))

pn_movement_visualization + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)



## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_pn_movement.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
embryo_aspect_ratio = embryo_scores[, c(2,3,9)]
embryo_aspect_ratio = drop_na(embryo_aspect_ratio)
embryo_aspect_ratio$embryo_aspect_ratio = as.numeric(embryo_aspect_ratio$embryo_aspect_ratio)


## -----------------------------------------------
emb_aspect_ratio_sample_sizes = calculate_sample_sizes(embryo_scores, "embryo_aspect_ratio")
sample_size_Nlab = sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(emb_aspect_ratio_sample_sizes$c48_hybrid,
                           emb_aspect_ratio_sample_sizes$cel_hybrid,
                           emb_aspect_ratio_sample_sizes$cre_hybrid,
                           emb_aspect_ratio_sample_sizes$csp48,
                          emb_aspect_ratio_sample_sizes$em464,
                          emb_aspect_ratio_sample_sizes$jk574, emb_aspect_ratio_sample_sizes$lkc28),")",sep="")

emb_aspect_ratio_descriptors = list("embryo_type")
create_factors(embryo_aspect_ratio, emb_aspect_ratio_descriptors)


## -----------------------------------------------
embryo_aspect_ratio_visualization = ggplot(embryo_aspect_ratio,
                                mapping = aes(x = embryo_type, y = embryo_aspect_ratio))

embryo_aspect_ratio_visualization + geom_boxplot()+ geom_jitter(width = 0.2) + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)



## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_embryo_aspect_ratio_jitter.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
centrosome_detachment = embryo_scores[, c(2,3,10)]
centrosome_detachment = drop_na(centrosome_detachment)
centrosome_detachment_grouped =  centrosome_detachment %>%
  group_by(embryo_type, centrosomal_detachment) %>% summarise(Count = n())


## -----------------------------------------------
centrosome_detachment_sample_sizes = calculate_sample_sizes(centrosome_detachment, "centrosomal_detachment")
sample_size_Nlab = sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(polar_body_sample_sizes$c48_hybrid,
                           polar_body_sample_sizes$cel_hybrid,
                           polar_body_sample_sizes$cre_hybrid,
                           polar_body_sample_sizes$csp48,
                          polar_body_sample_sizes$em464,
                          polar_body_sample_sizes$jk574, polar_body_sample_sizes$lkc28),")",sep="")

centrosome_detachment_descriptors = list("embryo_type")
create_factors(centrosome_detachment_grouped, centrosome_detachment_descriptors)


## -----------------------------------------------
centrosome_detachment_visualization = ggplot(centrosome_detachment_grouped,
                                mapping = aes(x = embryo_type, y = Count, fill = centrosomal_detachment))

centrosome_detachment_visualization + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)



## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_centrosomal_detachment.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
centrosome_positioning = embryo_scores[, c(2,3,11)]
centrosome_positioning = drop_na(centrosome_positioning)
centrosome_positioning_grouped =  centrosome_positioning %>%
  group_by(embryo_type, centrosome.positioning.defect.at.PNM) %>% summarise(Count = n())


## -----------------------------------------------
centrosome_positioning_sample_sizes = calculate_sample_sizes(centrosome_positioning, "centrosome.positioning.defect.at.PNM")
sample_size_Nlab = sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(centrosome_positioning_sample_sizes$c48_hybrid,
                           centrosome_positioning_sample_sizes$cel_hybrid,
                           centrosome_positioning_sample_sizes$cre_hybrid,
                           centrosome_positioning_sample_sizes$csp48,
                          centrosome_positioning_sample_sizes$em464,
                          centrosome_positioning_sample_sizes$jk574, centrosome_positioning_sample_sizes$lkc28),")",sep="")

centrosome_positioning_descriptors = list("embryo_type")
create_factors(centrosome_positioning_grouped, centrosome_positioning_descriptors)


## -----------------------------------------------
centrosome_positioning_visualization = ggplot(centrosome_positioning_grouped,
                                mapping = aes(x = embryo_type, y = Count, fill = centrosome.positioning.defect.at.PNM))

centrosome_positioning_visualization + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)



## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_centrosomal_positioning_at_PN.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
extra_cytoplasmic_clearings = embryo_scores[, c(2,3,12)]
extra_cytoplasmic_clearings = drop_na(extra_cytoplasmic_clearings)
extra_cytoplasmic_grouped =  extra_cytoplasmic_clearings %>%
  group_by(embryo_type, extra_cytoplasmic_clearings) %>% summarise(Count = n())


## -----------------------------------------------
extra_clearing_sample_sizes = calculate_sample_sizes(extra_cytoplasmic_clearings, "extra_cytoplasmic_clearings")
sample_size_Nlab = sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(polar_body_sample_sizes$c48_hybrid,
                           polar_body_sample_sizes$cel_hybrid,
                           polar_body_sample_sizes$cre_hybrid,
                           polar_body_sample_sizes$csp48,
                          polar_body_sample_sizes$em464,
                          polar_body_sample_sizes$jk574, polar_body_sample_sizes$lkc28),")",sep="")

extra_cytoplasmic_clearings_descriptors = list("embryo_type")
create_factors(extra_cytoplasmic_grouped, extra_cytoplasmic_clearings_descriptors)


## -----------------------------------------------
extra_cytoplasmic_clearings_visualization = ggplot(extra_cytoplasmic_grouped,
                                mapping = aes(x = embryo_type, y = Count, fill = extra_cytoplasmic_clearings))

extra_cytoplasmic_clearings_visualization + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)



## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_centrosomal_positioning_at_PN.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
spindle_orientation = embryo_scores[, c(2,3,13)]
spindle_orientation = drop_na(spindle_orientation)
spindle_orientation_grouped =  spindle_orientation %>%
  group_by(embryo_type, spindle_orientation_defects) %>% summarise(Count = n())


## -----------------------------------------------
spindle_orientation_defects_sample_sizes = calculate_sample_sizes(spindle_orientation, "spindle_orientation_defects")
sample_size_Nlab = sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(spindle_orientation_defects_sample_sizes$c48_hybrid,
                           spindle_orientation_defects_sample_sizes$cel_hybrid,
                           spindle_orientation_defects_sample_sizes$cre_hybrid,
                           spindle_orientation_defects_sample_sizes$csp48,
                          spindle_orientation_defects_sample_sizes$em464,
                          spindle_orientation_defects_sample_sizes$jk574, spindle_orientation_defects_sample_sizes$lkc28),")",sep="")

spindle_orientation_descriptors = list("embryo_type")
create_factors(spindle_orientation_grouped, spindle_orientation_descriptors)


## -----------------------------------------------
spindle_orientation_visualization = ggplot(spindle_orientation_grouped,
                                mapping = aes(x = embryo_type, y = Count, fill = spindle_orientation_defects))

spindle_orientation_visualization + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)



## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_spindle_orientation.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
four_cell_cytokinesis = embryo_scores[, c(2,3,14)]
four_cell_cytokinesis = drop_na(four_cell_cytokinesis)
four_cell_cytokinesis_grouped =  four_cell_cytokinesis %>%
  group_by(embryo_type, cytokinesis_failure_4.cell_stage) %>% summarise(Count = n())


## -----------------------------------------------
four_cell_cytokinesis_sample_sizes = calculate_sample_sizes(four_cell_cytokinesis, "cytokinesis_failure_4.cell_stage")
sample_size_Nlab = sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(four_cell_cytokinesis_sample_sizes$c48_hybrid,
                           four_cell_cytokinesis_sample_sizes$cel_hybrid,
                           four_cell_cytokinesis_sample_sizes$cre_hybrid,
                           four_cell_cytokinesis_sample_sizes$csp48,
                          four_cell_cytokinesis_sample_sizes$em464,
                          four_cell_cytokinesis_sample_sizes$jk574, four_cell_cytokinesis_sample_sizes$lkc28),")",sep="")

four_cell_cytokinesis_descriptors = list("embryo_type")
create_factors(spindle_orientation_grouped, four_cell_cytokinesis_descriptors)


## -----------------------------------------------
four_cell_cytokinesis_visualization = ggplot(four_cell_cytokinesis_grouped,
                                mapping = aes(x = embryo_type, y = Count, fill = cytokinesis_failure_4.cell_stage))

four_cell_cytokinesis_visualization + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)



## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_four_cell_cytokinesis_failure.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
excess_contractility = embryo_scores[, c(2,3,15)]
excess_contractility = drop_na(excess_contractility)
excess_contractility_grouped =  excess_contractility %>%
  group_by(embryo_type, excess_contractility) %>% summarise(Count = n())


## -----------------------------------------------
excess_contractility_sample_sizes = calculate_sample_sizes(excess_contractility, "cytokinesis_failure_4.cell_stage")
sample_size_Nlab = sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(excess_contractility_sample_sizes$c48_hybrid,
                           excess_contractility_sample_sizes$cel_hybrid,
                           excess_contractility_sample_sizes$cre_hybrid,
                           excess_contractility_sample_sizes$csp48,
                          excess_contractility_sample_sizes$em464,
                          excess_contractility_sample_sizes$jk574, excess_contractility_sample_sizes$lkc28),")",sep="")

excess_contractility_descriptors = list("embryo_type")
create_factors(excess_contractility_grouped, excess_contractility_descriptors)


## -----------------------------------------------
excess_contractility_visualization = ggplot(excess_contractility_grouped,
                                mapping = aes(x = embryo_type, y = Count, fill = excess_contractility))

excess_contractility_visualization + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)



## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_excess_contractility.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
two_cell_junction = embryo_scores[, c(2,3,16)]
two_cell_junction = drop_na(two_cell_junction)
two_cell_junction_grouped =  two_cell_junction %>%
  group_by(embryo_type, rounded.2.cell.junction) %>% summarise(Count = n())


## -----------------------------------------------
two_cell_junction_sample_sizes = calculate_sample_sizes(two_cell_junction, "rounded.2.cell.junction")
sample_size_Nlab = sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(two_cell_junction_sample_sizes$c48_hybrid,
                           two_cell_junction_sample_sizes$cel_hybrid,
                           two_cell_junction_sample_sizes$cre_hybrid,
                           two_cell_junction_sample_sizes$csp48,
                          two_cell_junction_sample_sizes$em464,
                          two_cell_junction_sample_sizes$jk574, two_cell_junction_sample_sizes$lkc28),")",sep="")

two_cell_junction_descriptors = list("embryo_type")
create_factors(two_cell_junction_grouped, two_cell_junction_descriptors)


## -----------------------------------------------
two_cell_junction_visualization = ggplot(two_cell_junction_grouped,
                                mapping = aes(x = embryo_type, y = Count, fill = rounded.2.cell.junction))

two_cell_junction_visualization + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)



## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_two_cell_junction.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
two_cell_evenness = embryo_scores[, c(2,3,17)]
two_cell_evenness = drop_na(two_cell_evenness)
two_cell_evenness_grouped =  two_cell_evenness %>%
  group_by(embryo_type, even_two_cells) %>% summarise(Count = n())


## -----------------------------------------------
two_cell_evenness_sample_sizes = calculate_sample_sizes(two_cell_junction, "rounded.2.cell.junction")
sample_size_Nlab = sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(two_cell_evenness_sample_sizes$c48_hybrid,
                           two_cell_evenness_sample_sizes$cel_hybrid,
                           two_cell_evenness_sample_sizes$cre_hybrid,
                           two_cell_evenness_sample_sizes$csp48,
                          two_cell_evenness_sample_sizes$em464,
                          two_cell_evenness_sample_sizes$jk574, two_cell_evenness_sample_sizes$lkc28),")",sep="")

two_cell_evenness_descriptors = list("embryo_type")
create_factors(two_cell_evenness_grouped, two_cell_evenness_descriptors)


## -----------------------------------------------
two_cell_evenness_visualization = ggplot(two_cell_evenness_grouped,
                                mapping = aes(x = embryo_type, y = Count, fill = even_two_cells))

two_cell_evenness_visualization + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)



## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_two_cell_evenness.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
split_two_cell_nucleus = embryo_scores[, c(2,3,18)]
split_two_cell_nucleus = drop_na(split_two_cell_nucleus)
split_two_cell_grouped =  split_two_cell_nucleus %>%
  group_by(embryo_type, split_nucleus_two_cell) %>% summarise(Count = n())


## -----------------------------------------------
split_two_cell_nuc_sample_sizes = calculate_sample_sizes(split_two_cell_nucleus, "split_nucleus_two_cell")
sample_size_Nlab = sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(split_two_cell_nuc_sample_sizes$c48_hybrid,
                           split_two_cell_nuc_sample_sizes$cel_hybrid,
                           split_two_cell_nuc_sample_sizes$cre_hybrid,
                           split_two_cell_nuc_sample_sizes$csp48,
                          split_two_cell_nuc_sample_sizes$em464,
                          split_two_cell_nuc_sample_sizes$jk574, split_two_cell_nuc_sample_sizes$lkc28),")",sep="")

split_two_cell_nucleus_descriptors = list("embryo_type")
create_factors(split_two_cell_grouped, split_two_cell_nucleus_descriptors)


## -----------------------------------------------
split_two_cell_nucleus_visualization = ggplot(split_two_cell_grouped,
                                mapping = aes(x = embryo_type, y = Count, fill = split_nucleus_two_cell))

split_two_cell_nucleus_visualization + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)



## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_split_nucleus_two_cell.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
cytokinesis_two_cell = embryo_scores[, c(2,3,19)]
cytokinesis_two_cell = drop_na(cytokinesis_two_cell)
cytokinesis_two_cell_grouped =  cytokinesis_two_cell %>%
  group_by(embryo_type, cytokinesis_failure_two.cell) %>% summarise(Count = n())


## -----------------------------------------------
cytokinesis_two_cell_failure_sample_sizes = calculate_sample_sizes(cytokinesis_two_cell, "cytokinesis_failure_two.cell")
sample_size_Nlab = sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(cytokinesis_two_cell_failure_sample_sizes$c48_hybrid,
                           cytokinesis_two_cell_failure_sample_sizes$cel_hybrid,
                           cytokinesis_two_cell_failure_sample_sizes$cre_hybrid,
                           cytokinesis_two_cell_failure_sample_sizes$csp48,
                          cytokinesis_two_cell_failure_sample_sizes$em464,
                          cytokinesis_two_cell_failure_sample_sizes$jk574, cytokinesis_two_cell_failure_sample_sizes$lkc28),")",sep="")

cytokinesis_two_cell_descriptors = list("embryo_type")
create_factors(cytokinesis_two_cell_grouped, cytokinesis_two_cell_descriptors)


## -----------------------------------------------
cytokinesis_two_cell_visualization = ggplot(cytokinesis_two_cell_grouped,
                                mapping = aes(x = embryo_type, y = Count, fill = cytokinesis_failure_two.cell))

cytokinesis_two_cell_visualization + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)



## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_cytokinesis_two_cell.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
four_cell_crosseyed = embryo_scores[, c(2,3,20)]
four_cell_crosseyed = drop_na(four_cell_crosseyed)
four_cell_crosseyed_grouped =  four_cell_crosseyed %>%
  group_by(embryo_type,four_cells_cross.eyes) %>% summarise(Count = n())


## -----------------------------------------------
four_cell_crosseyed_sample_sizes = calculate_sample_sizes(four_cell_crosseyed, "four_cells_cross.eyes")
sample_size_Nlab = sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(four_cell_crosseyed_sample_sizes$c48_hybrid,
                           four_cell_crosseyed_sample_sizes$cel_hybrid,
                           four_cell_crosseyed_sample_sizes$cre_hybrid,
                           four_cell_crosseyed_sample_sizes$csp48,
                          four_cell_crosseyed_sample_sizes$em464,
                          four_cell_crosseyed_sample_sizes$jk574, four_cell_crosseyed_sample_sizes$lkc28),")",sep="")

four_cell_crosseyes_descriptors = list("embryo_type")
create_factors(four_cell_crosseyed_grouped, four_cell_crosseyes_descriptors)


## -----------------------------------------------
four_cell_crosseyed_visualization = ggplot(four_cell_crosseyed_grouped,
                                mapping = aes(x = embryo_type, y = Count, fill = four_cells_cross.eyes))

four_cell_crosseyed_visualization + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)



## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_four_cell_crosseyes.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
four_cell_polarity = embryo_scores[, c(2,3,21)]
four_cell_polarity = drop_na(four_cell_polarity)
four_cell_polarity_grouped =  four_cell_polarity %>%
  group_by(embryo_type,four_cell_polarity_defect) %>% summarise(Count = n())


## -----------------------------------------------
four_cell_polarity_sample_sizes = calculate_sample_sizes(four_cell_polarity, "four_cell_polarity_defect")
sample_size_Nlab = sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(four_cell_polarity_sample_sizes$c48_hybrid,
                           four_cell_polarity_sample_sizes$cel_hybrid,
                           four_cell_polarity_sample_sizes$cre_hybrid,
                           four_cell_polarity_sample_sizes$csp48,
                          four_cell_polarity_sample_sizes$em464,
                          four_cell_polarity_sample_sizes$jk574, four_cell_polarity_sample_sizes$lkc28),")",sep="")

four_cell_polarity_descriptors = list("embryo_type")
create_factors(four_cell_crosseyed_grouped, four_cell_polarity_descriptors)


## -----------------------------------------------
four_cell_polarity_visualization = ggplot(four_cell_polarity_grouped,
                                mapping = aes(x = embryo_type, y = Count, fill = four_cell_polarity_defect))

four_cell_polarity_visualization + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)



## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_four_cell_polarity.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
par_two_spindle_orientation = embryo_scores[, c(2,3,22)]
par_two_spindle_orientation = drop_na(par_two_spindle_orientation)
par_two_spindle_grouped =  par_two_spindle_orientation %>%
  group_by(embryo_type,par.2.like.spindle.orientation) %>% summarise(Count = n())


## -----------------------------------------------
par_two_spindle_sample_sizes = calculate_sample_sizes(par_two_spindle_orientation, "par.2.like.spindle.orientation")
sample_size_Nlab = sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(par_two_spindle_sample_sizes$c48_hybrid,
                           par_two_spindle_sample_sizes$cel_hybrid,
                           par_two_spindle_sample_sizes$cre_hybrid,
                           par_two_spindle_sample_sizes$csp48,
                          par_two_spindle_sample_sizes$em464,
                          par_two_spindle_sample_sizes$jk574, par_two_spindle_sample_sizes$lkc28),")",sep="")

par_two_spindle_descriptors = list("embryo_type")
create_factors(par_two_spindle_grouped, par_two_spindle_descriptors)


## -----------------------------------------------
par_two_spindle_visualization = ggplot(par_two_spindle_grouped,
                                mapping = aes(x = embryo_type, y = Count, fill = par.2.like.spindle.orientation))

par_two_spindle_visualization + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)



## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_par-2_spindle_orientation.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
par_six_spindle_orientation = embryo_scores[, c(2,3,23)]
par_six_spindle_orientation = drop_na(par_six_spindle_orientation)
par_six_spindle_grouped =  par_six_spindle_orientation %>%
  group_by(embryo_type,par.6.like.spindle.orientation) %>% summarise(Count = n())


## -----------------------------------------------
par_six_spindle_sample_sizes = calculate_sample_sizes(par_six_spindle_orientation, "par.6.like.spindle.orientation")
sample_size_Nlab = sample_size_Nlab = paste(c("c48_hybrid",
                                "cel_hybrid", "cre_hybrid", "csp48", "crem", "cel", "cbn"),"\n(N=",
                         c(par_six_spindle_sample_sizes$c48_hybrid,
                           par_six_spindle_sample_sizes$cel_hybrid,
                           par_six_spindle_sample_sizes$cre_hybrid,
                           par_six_spindle_sample_sizes$csp48,
                          par_six_spindle_sample_sizes$em464,
                          par_six_spindle_sample_sizes$jk574, par_six_spindle_sample_sizes$lkc28),")",sep="")

par_six_spindle_descriptors = list("embryo_type")
create_factors(par_six_spindle_grouped, par_six_spindle_descriptors)


## -----------------------------------------------
par_six_spindle_visualization = ggplot(par_six_spindle_grouped,
                                mapping = aes(x = embryo_type, y = Count, fill = par.6.like.spindle.orientation))

par_six_spindle_visualization + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= sample_size_Nlab)



## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_par-6_spindle_orientation.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
evo_adjusted_phenotypes = read_csv(paste0(blinded_embryo_scores_file_directory,blinded_embryo_scores_for_heat_map))

## -----------------------------------------------
embryo_scores = pivot_longer(data = evo_adjusted_phenotypes,
                             cols = "c48_hybrid":"lkc28",
                             names_to = "embryo_type", 
                             values_to = "phenotype_score")

embryo_phenotypes = list(unique(embryo_scores$phenotype_name))
embryo_types = list(unique(embryo_scores$embryo_type))

general_embryo_types = list()
for (embryo in embryo_scores$embryo_type){
  if (grepl("hybrid", embryo) == TRUE){
    embryo_type = "hybrid"
  }else{
    embryo_type = "wild-type"
  }
  general_embryo_types = append(general_embryo_types, embryo_type)
}

embryo_scores$general_embryo_type = general_embryo_types


## -----------------------------------------------
embryo_scores$general_embryo_type <- factor(embryo_scores$general_embryo_type,
          levels = unique(embryo_scores$general_embryo_type))


## -----------------------------------------------
cell_division_ordered_phenotypes = c("irregular_polarbody", "uneven_sized_nuclei_pronuclear-meeting", "pronuclei_fail_to_meet",
                                        "misshappen_nuclei", "abnormal_pronuclei_movement","centrosomal_detachment",
                                        "centrosome_positioning_defect_at_PNM", "extra_cytoplasmic_clearings",
                                        "spindle_orientation_defects", "binary_embryo_aspect_ratio", "excess_contractility",
                                        "rounded_2-cell_junction", "even_two_cells",
                                        "split_nucleus_two_cell", "cytokinesis_failure_two-cell", "four_cells_cross-eyes",
                                        "four_cell_polarity_defect","par-2-like-spindle-orientation","par-6-like-spindle-orientation",
                                     "cytokinesis_failure_4 cell_stage") 
embryo_ordered_by_phylogeny = c("lkc28", "csp48", "em464", "jk574", "c48_hybrid", "cre_hybrid", "cel_hybrid")


## -----------------------------------------------
embryo_scores_phenotypes_removed = subset(embryo_scores, phenotype_name!= "pronuclei_fail_to_meet" & phenotype_name != "abnormal_pronuclei_movement" & phenotype_name != "extra_cytoplasmic_clearings" & phenotype_name != "excess_contractility" & phenotype_name != "cytokinesis_failure_two-cell" & phenotype_name != "par-6-like-spindle-orientation")

ordered_phenotypes_winnowed = c("irregular_polarbody", "uneven_sized_nuclei_pronuclear-meeting",
                                        "misshappen_nuclei","centrosomal_detachment",
                                        "centrosome_positioning_defect_at_PNM",
                                        "spindle_orientation_defects", "binary_embryo_aspect_ratio",
                                        "rounded_2-cell_junction", "even_two_cells",
                                        "split_nucleus_two_cell", "four_cells_cross-eyes",
                                        "four_cell_polarity_defect","par-2-like-spindle-orientation",
                                     "cytokinesis_failure_4 cell_stage") 


## -----------------------------------------------
wt_phenotype_scores = subset(embryo_scores, grepl("hybrid", embryo_type) == FALSE)
hybrid_phenotype_scores = subset(embryo_scores, grepl("hybrid", embryo_type) == TRUE)

wt_phenotype_scores_winnowed = subset(embryo_scores_phenotypes_removed, grepl("hybrid", embryo_type) == FALSE)
hybrid_phenotype_scores_winnowed = subset(embryo_scores_phenotypes_removed, grepl("hybrid", embryo_type) == TRUE)



## -----------------------------------------------
phenotype_score_range = range(c((wt_phenotype_scores["phenotype_score"]), (hybrid_phenotype_scores["phenotype_score"])))

phenotype_score_range_winnowed = range(c((wt_phenotype_scores_winnowed["phenotype_score"]), (hybrid_phenotype_scores_winnowed["phenotype_score"])))


## -----------------------------------------------
wt_phenotype_scores_heatmap = ggplot(data = wt_phenotype_scores, mapping = aes(x = embryo_type, y = phenotype_name, fill = phenotype_score)) 
wt_phenotype_scores_heatmap + geom_tile() + theme_cowplot(14) + scale_fill_gradient2(low = "white", mid = "deepskyblue", high = "darkblue",midpoint = mean(phenotype_score_range), breaks = seq(0,1,0.25), limits=c(floor(phenotype_score_range[1]), ceiling(phenotype_score_range[2])))


## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_wt_crosses_heatmap_no_aspect_ratio.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
hybrid_phenotype_scores_heatmap = ggplot(data = hybrid_phenotype_scores, mapping = aes(x = embryo_type, y = phenotype_name, fill = phenotype_score)) 
hybrid_phenotype_scores_heatmap + geom_tile() + theme_cowplot(14) + scale_fill_gradient2(low = "white", mid = "deepskyblue", high = "darkblue",midpoint = mean(phenotype_score_range), breaks = seq(0,1,0.25), limits=c(floor(phenotype_score_range[1]), ceiling(phenotype_score_range[2])))

## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_hybrid_crosses_heatmap_no_aspect_ratio.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
phenotype_scores_heatmap = ggplot(data = embryo_scores, mapping = aes(x = embryo_type, y = phenotype_name, fill = phenotype_score)) 
phenotype_scores_heatmap + geom_tile() + theme_cowplot(14) + scale_fill_gradient2(low = "snow", mid = "deepskyblue", high = "darkblue", midpoint = mean(phenotype_score_range), breaks = seq(0,1,0.25), limits=c(floor(phenotype_score_range[1]), ceiling(phenotype_score_range[2]))) + scale_y_discrete(limits = cell_division_ordered_phenotypes) + scale_x_discrete(limits = embryo_ordered_by_phylogeny)


## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_hybrid_crosses_heatmap.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## -----------------------------------------------
phenotype_scores_heatmap = ggplot(data = embryo_scores_phenotypes_removed, mapping = aes(x = embryo_type, y = phenotype_name, fill = phenotype_score)) 
phenotype_scores_heatmap + geom_tile() + theme_cowplot(14) + scale_fill_gradient2(low = "snow", mid = "deepskyblue", high = "darkblue", midpoint = mean(phenotype_score_range_winnowed), breaks = seq(0,1,0.25), limits=c(floor(phenotype_score_range_winnowed[1]), ceiling(phenotype_score_range_winnowed[2]))) + scale_y_discrete(limits = ordered_phenotypes_winnowed) + scale_x_discrete(limits = embryo_ordered_by_phylogeny)


## -----------------------------------------------
ggsave(paste0(blinded_embryo_scores_file_directory,"evo_crosses_phenotypes_winnowed.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)

