# Clustered heat map of sequence divergence scores

# much of this code was based off this instructional guide: https://jcoliver.github.io/learn-r/008-ggplot-dendrograms-and-heatmaps.html

library(ggplot2); library(dplyr); library(reshape); library(Hmisc);
library(tidyverse); library(stats); library(EnvStats); library(cowplot); library(devtools);
library(wbData); library(stringr); library(ggdendro); library(grid); library(data.table)

get_file_info = function(trailingOnly = TRUE, asValues = TRUE) {
  data_file_info = commandArgs(TRUE)
  if (length(data_file_info) < 4){
    stop("Provide the file path to functions,
         data file dir, data file name, simple mine file name,
         heatmap export name, eg. dev/ sequence_divergence_scores_data_dir/ sequence_div_data.ext simple_mine_data.ext seq_evo_divergence_heat_map_from_python.ext")
  } else {
    return(data_file_info)
  }
}

#Get file info
sequence_divergence_file_info = get_file_info()
function_directory = sequence_divergence_file_info[1]
sequence_divergence_file_directory = sequence_divergence_file_info[2]
sequence_divergence_filename = sequence_divergence_file_info[3]
simple_mine_data = sequence_divergence_file_info[4]
sequence_divergence_heatmap = sequence_divergence_file_info[5]

## Import global user functions
setwd(function_directory)
source("get_functions.R")
user_functions = list("select_variables_reformat_long.R","create_factors_sample_sizes.R",
                      "centriole_localization_data_handling.R")

get_functions(function_directory, user_functions)

## import divergence data

cbn_sequence_divergence_scores = read.csv(paste0(sequence_divergence_file_directory,
                                                 sequence_divergence_filename))
## remove cbn scores

cbn_sequence_divergence_scores = cbn_sequence_divergence_scores %>% slice(2:5)

## get geneids for all transcript names
## need to create simplemine file with transcript names in the format "CELE_transcript"
## Collect the following data: wormbase geneID, public name, sequence name, uniProt,
## Gene ontology association, Concise description

simple_mine_data = read.csv(paste0(sequence_divergence_file_directory, simple_mine_data))

gene_public_names = simple_mine_data$public_name
gene_public_names_with_title = c("species_name", gene_public_names)
colnames(cbn_sequence_divergence_scores) = gene_public_names_with_title

## create clustered heat map
species_ordered_by_relatedness = c("CSP48", "CSINI", "CREMA", "CELEG")

cbn_divergence_scores_transposed = transpose(cbn_sequence_divergence_scores)
rownames(cbn_divergence_scores_transposed) = colnames(cbn_sequence_divergence_scores)
colnames(cbn_divergence_scores_transposed) = cbn_sequence_divergence_scores$species_name
cbn_divergence_scores_transposed = cbn_divergence_scores_transposed[-c(1),]
seq_evo_divergence_scores = as.matrix(cbn_divergence_scores_transposed)
seq_evo_divergence_matrix = as.dendrogram(hclust(d = dist(x = seq_evo_divergence_scores)))

dendro_plot = ggdendrogram(data = seq_evo_divergence_matrix, rotate = TRUE)
dendro_plot + theme(axis.text.y = element_text(size = 6))

ggsave(paste0(sequence_divergence_file_directory,"cbn_sequence_divergence_dendrogram.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)

## create heat map
gene_names = rownames(cbn_divergence_scores_transposed)
cbn_divergence_scores_transposed$gene_names = gene_names
cbn_divergence_scores_for_heatmap = pivot_longer(cbn_divergence_scores_transposed,
                                                 cols = -c(gene_names),
                                                 names_to = "species",
                                                 values_to = "divergence_from_cbn")
cbn_divergence_scores_for_heatmap$divergence_from_cbn = as.numeric(cbn_divergence_scores_for_heatmap$divergence_from_cbn)
cbn_divergence_heatmap = ggplot(cbn_divergence_scores_for_heatmap, aes(x = species, y = gene_names, fill = divergence_from_cbn))
cbn_divergence_heatmap + geom_tile((aes(fill = divergence_from_cbn))) + theme(axis.text.y = element_text(size = 6)) + scale_x_discrete(limits = species_ordered_by_relatedness)

## combine the dendrogram and the heatmap

gene_order = order.dendrogram(seq_evo_divergence_matrix)
cbn_divergence_scores_for_heatmap$gene_names <- factor(cbn_divergence_scores_for_heatmap$gene_names,
                                                      levels = cbn_divergence_scores_transposed$gene_names[gene_order],
                                                      ordered = TRUE)

cbn_divergence_heatmap = ggplot(cbn_divergence_scores_for_heatmap,
                                aes(x = species, y = gene_names,
                                    fill = divergence_from_cbn)) + geom_tile((aes(fill = divergence_from_cbn))) + theme(axis.text.y = element_text(size = 6),
                                                                                                                        legend.position = "top") + theme_bw() + scale_x_discrete(limits = species_ordered_by_relatedness) + scale_fill_gradient(low = "snow", high = "darkorchid")
ggsave(paste0(sequence_divergence_file_directory,"cbn_sequence_divergence_heatmap.pdf"), width = 20,
       height = 20, units = "cm",
       device = "pdf", dpi = 300)

