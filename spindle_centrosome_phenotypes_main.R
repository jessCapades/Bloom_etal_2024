## ----setup, include=FALSE---------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---------------------------------------------------------------------------------
library(ggplot2); library(dplyr); library(reshape); library(Hmisc);
library(tidyverse); library(stats); library(EnvStats); library(cowplot)

## ---------------------------------------------------------------------------------
get_file_info = function(trailingOnly = TRUE, asValues = TRUE) {
  data_file_info = commandArgs(TRUE)
  if (length(data_file_info) < 3){
    stop("Provide the file path to functions,
         data file dir, and filename, eg. dev/ spindle_centrosome_data_dir/ spindle_centrosome_data.ext")
  } else {
    return(data_file_info)
  }
}

spindle_centrosome_file_info = get_file_info()
function_directory = spindle_centrosome_file_info[1]
spindle_centrosome_file_directory = spindle_centrosome_file_info[2]
spindle_centrosome_filename = spindle_centrosome_file_info[3]
## ---------------------------------------------------------------------------------
spindle_phenotypes = read.csv(paste0(spindle_centrosome_file_directory, spindle_centrosome_filename))
## Import global user functions
setwd(function_directory)
source("get_functions.R")

user_functions = list("select_variables_reformat_long.R","create_factors_sample_sizes.R",
                      "centriole_localization_data_handling.R")

get_functions(function_directory, user_functions)


## ---------------------------------------------------------------------------------
spindle_phenotypes = spindle_phenotypes[, c(1:2, 8:10)]
spindle_phenotypes = subset(spindle_phenotypes, embryo_id !="", select = c(1:5))
spindle_phenotypes = spindle_phenotypes[!grepl("eepi11324_od3701", spindle_phenotypes$embryo_id), ]

centriole_phenotypes = spindle_phenotypes[, c(1:3)]
centriole_phenotypes = drop_na(centriole_phenotypes)
spindle_broad_phenotypes = spindle_phenotypes[, c(1:2,4)]
spindle_broad_phenotypes = drop_na(spindle_broad_phenotypes)
meiotic_spindle_capture = spindle_phenotypes[, c(1:2, 5)]
meiotic_spindle_capture = drop_na(meiotic_spindle_capture)

spindle_phenotype_descriptors = c("embryo_type")

phenotype_descriptions = list(centriole_phenotypes, spindle_broad_phenotypes, meiotic_spindle_capture)

for(description in phenotype_descriptions) {
  create_factors(description, spindle_phenotype_descriptors)
}

centriole_phenotypes_grouped = centriole_phenotypes %>%
  group_by(embryo_type, early_centriole_separation) %>%
  summarise(Count = n())

spindle_broad_phenotypes_grouped = spindle_broad_phenotypes %>%
  group_by(embryo_type, misformed_spindle) %>%
  summarise(Count = n())

meiotic_spindle_grouped = meiotic_spindle_capture %>%
  group_by(embryo_type, meiotic_spindle_capture) %>%
  summarise(Count = n())


## ---------------------------------------------------------------------------------
centriole_detachment_sample_sizes = calculate_sample_sizes(centriole_phenotypes, "early_centriole_separation")

centriole_detachment_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                            Hybrid"),"\n(N=", c(centriole_detachment_sample_sizes$od3701,
                                                centriole_detachment_sample_sizes$lkc28,
                                          centriole_detachment_sample_sizes$hybrid),")",sep="")

centriole_position_vis = ggplot(centriole_phenotypes_grouped, mapping = aes(embryo_type, Count, fill = early_centriole_separation))
centriole_position_vis + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels= centriole_detachment_Nlab)

ggsave(paste0(spindle_centrosome_file_directory,"sirTub_centriole_positioning.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
spindle_formation_sample_sizes = calculate_sample_sizes(spindle_broad_phenotypes, "misformed_spindle")

spindle_formation_sample_sizes_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                            Hybrid"),"\n(N=", c(spindle_formation_sample_sizes$od3701,
                                                spindle_formation_sample_sizes$lkc28,
                                          spindle_formation_sample_sizes$hybrid),")",sep="")
spindle_formation_vis = ggplot(spindle_broad_phenotypes_grouped, mapping = aes(embryo_type, Count, fill = misformed_spindle))

spindle_formation_vis + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels = spindle_formation_sample_sizes_Nlab)

ggsave(paste0(spindle_centrosome_file_directory,"sirTub_spindle_broad_phenotypes.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
meiotic_capture_sample_sizes = calculate_sample_sizes(meiotic_spindle_capture, "meiotic_spindle_capture")

meiotic_capture_sample_sizes_Nlab = paste(c("C. elegans Wild-type",
                                "C. brenneri Wild-type",
                                "C. brenneri(f) x C. elegans(m) 
                            Hybrid"),"\n(N=", c(meiotic_capture_sample_sizes$od3701,
                                                meiotic_capture_sample_sizes$lkc28,
                                          meiotic_capture_sample_sizes$hybrid),")",sep="")

meiotic_capture_vis = ggplot(meiotic_spindle_grouped, mapping = aes(embryo_type, Count, fill = meiotic_spindle_capture))

meiotic_capture_vis + geom_bar(stat = "identity", position = "fill") + theme_cowplot(14) + scale_x_discrete(labels = meiotic_capture_sample_sizes_Nlab)

ggsave(paste0(spindle_centrosome_file_directory,"sirTub_meiotic_spindle_capture.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)

