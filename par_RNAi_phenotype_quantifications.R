# par RNAi phenotype quantifications

#Import libraries
library(ggplot2); library(dplyr); library(reshape); library(Hmisc); library(tidyverse); library(stats);
library(EnvStats); library(cowplot)
## Import global user functions

rm(list=setdiff(ls(), "par_phenotype_measurements"))

get_file_info = function(trailingOnly = TRUE, asValues = TRUE) {
  data_file_info = commandArgs(TRUE)
  if (length(data_file_info) < 4){
    stop("Provide the file path to functions,
         data file dir, and filenames, eg. dev/ par_phenotypes_data_dir/ cbn_par_phenotypes_data.ext cel_par_phenotypes_data.ext")
  } else {
    return(data_file_info)
  }
}

par_phenotypes_file_info = get_file_info()
function_directory = par_phenotypes_file_info[1]
par_phenotypes_file_directory = par_phenotypes_file_info[2]
cbrenneri_phenotypes_filename = par_phenotypes_file_info[3]
celegans_phenotypes_filename = par_phenotypes_file_info[4]

setwd(function_directory)
source(paste0(function_directory,"select_variables_reformat_long.R"))
source(paste0(function_directory,"create_factors_sample_sizes.R"))
source(paste0(function_directory,"two_cell_size_comparisons.R"))
source(paste0(function_directory,"four_cell_phenotype_frequencies.R"))

par_phenotype_measurements = read.csv(paste0(par_phenotypes_file_directory, cbrenneri_phenotypes_filename))


par_rnai_phenotypes = par_phenotype_measurements %>%
  group_by(rnai, four_cell_phenotype, imaging_time)%>%
  summarise(Count = n())

data_descriptors = list("four_cell_phenotype", "imaging_time")
create_factors(par_rnai_phenotypes, data_descriptors = data_descriptors )

par_rnai_vis = ggplot(par_rnai_phenotypes, mapping = aes(x = imaging_time, y = Count, fill = four_cell_phenotype))

par_rnai_vis + geom_bar(position = "fill", stat = "identity") + theme_cowplot(14) + facet_wrap(~rnai)

ggsave(paste0(par_phenotypes_file_directory,"par_rnai_phenotypes.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)

cel_par_phenotype_measurements = read.csv(paste0(par_phenotypes_file_directory,celegans_phenotypes_filename))

cel_par_rnai_phenotypes = cel_par_phenotype_measurements %>%
  group_by(rnai_species, phenotype)%>%
  summarise(Count = n())

data_descriptors = list("phenotype")
create_factors(cel_par_rnai_phenotypes, data_descriptors = data_descriptors )

cel_par_rnai_vis = ggplot(cel_par_rnai_phenotypes, mapping = aes(x = rnai_species,
                                                                 y = Count, fill = phenotype))

cel_par_rnai_vis + geom_bar(position = "fill", stat = "identity") + theme_cowplot(14)

ggsave(paste0(par_phenotypes_file_directory,"cel_par_rnai_phenotypes.pdf"),
       width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)

