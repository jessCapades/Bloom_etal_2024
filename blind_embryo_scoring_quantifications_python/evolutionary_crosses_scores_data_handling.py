## pseudo code
## calculate normalized phenotype score
## combine to one dataframe

import pandas as pd
import numpy as np
from statistics import stdev
from math import isnan

# create new embryo_aspect_ratio_variable for later calculations

def calculate_wild_type_aspect_95_conf_int(phenotype_scores):
    
    embryo_aspect_ratios = phenotype_scores[['embryo_type','embryo_aspect_ratio']]
    #filter for scores of only wild-type embryos
    wild_type_embryos = ['lkc28', 'csp48', 'em464', 'jk574']
    
    wild_type_embryo_mask = embryo_aspect_ratios['embryo_type'].isin(wild_type_embryos)
    wild_type_embryo_aspect_ratios = embryo_aspect_ratios[wild_type_embryo_mask]
    wild_type_embryo_aspect_ratios_average = sum(wild_type_embryo_aspect_ratios['embryo_aspect_ratio'])/len(wild_type_embryo_aspect_ratios['embryo_aspect_ratio'])
    wild_type_embryo_aspect_ratios_standard_dev = stdev(wild_type_embryo_aspect_ratios['embryo_aspect_ratio'])
    
    #calculate data within 2SD range of wild-type average embryo aspect ratios
    wild_type_upper_interval = wild_type_embryo_aspect_ratios_average + 2*wild_type_embryo_aspect_ratios_standard_dev
    wild_type_lower_interval = wild_type_embryo_aspect_ratios_average - 2*wild_type_embryo_aspect_ratios_standard_dev
    wild_type_embryo_aspect_ratio_95confint = [wild_type_lower_interval, wild_type_upper_interval]
    
    return(wild_type_embryo_aspect_ratio_95confint)

def calculate_aspect_ratio_score(phenotype_scores):

    binary_aspect_ratio_scores = {}
    embryo_aspect_ratios = phenotype_scores[['embryo_aspect_ratio']]
    aspect_ratios = embryo_aspect_ratios.to_dict()
    wild_type_95_conf_int = calculate_wild_type_aspect_95_conf_int(phenotype_scores)
    
    binary_score = []
    indexed_aspect_ratios = aspect_ratios.get('embryo_aspect_ratio')

    for index in indexed_aspect_ratios.keys():

        aspect_ratio_score = indexed_aspect_ratios.get(index)

        if isnan(aspect_ratio_score):
            binary_score = 'empty'
        else:
            if wild_type_95_conf_int[0] < aspect_ratio_score < wild_type_95_conf_int[1]:
                binary_score = 0
            else:
                binary_score = 1
        binary_aspect_ratio_scores[index] = binary_score

    return(binary_aspect_ratio_scores)

## order binary embryo scores to match original dataframe
def append_binary_aspect_ratios_to_phenotype_scores (phenotype_scores, binary_aspect_ratio_scores):

    binary_aspect_ratios = binary_aspect_ratio_scores.items()
    sorted_binary_aspect_ratios = sorted(binary_aspect_ratios, key=lambda tup: tup[0])
    sorted_binary_aspect_ratio_scores = [x[1] for x in sorted_binary_aspect_ratios]
    phenotype_scores_with_binary_aspect_ratios = phenotype_scores.assign(binary_embryo_aspect_ratio = sorted_binary_aspect_ratio_scores)
    phenotype_scores_with_binary_aspect_ratios.replace('empty', np.nan, inplace=True)
    
    return(phenotype_scores_with_binary_aspect_ratios)

## calculate grouped phenotype scores and sample sizes
def calculate_phenotype_scores(phenotype_scores):
    
    phenotypes = phenotype_scores.columns.values.tolist()
    embryo_scores_per_phenotype = {}
    for phenotype in phenotypes[3:]:

        grouped_phenotype = phenotype_scores.groupby(['embryo_type'])[phenotype].sum()
        grouped_phenotype_lookup = grouped_phenotype.to_dict()
        embryo_scores_per_phenotype[phenotype] = grouped_phenotype_lookup

    return embryo_scores_per_phenotype

def calculate_sample_size(phenotype_scores):
    
    emb_phenotypes = phenotype_scores.columns.values.tolist()
    embryo_number_per_phenotype = {}
    for emb_phenotype in emb_phenotypes[3:]:
        embryo_numbers = phenotype_scores.groupby(['embryo_type'])[emb_phenotype].count()
        grouped_sample_size = embryo_numbers.to_dict()
        embryo_number_per_phenotype[emb_phenotype] = grouped_sample_size
    
    return embryo_number_per_phenotype

## calculate normalized phenotype scores

def calculate_phenotype_adjusted_score(embryo_scores_per_phenotype, embryo_sample_sizes):

    phenotype_adjusted_scores = {}
    for phenotype in embryo_scores_per_phenotype.keys():
        emb_phenotypes = embryo_scores_per_phenotype.get(phenotype)
        emb_sample_sizes = embryo_sample_sizes.get(phenotype)
        embryos_adjusted_scores = {}
        for embryo_type in emb_phenotypes.keys():
            species_score = emb_phenotypes.get(embryo_type)
            species_sample_size = emb_sample_sizes.get(embryo_type)
            species_adjusted_score = species_score/species_sample_size
            embryos_adjusted_scores[embryo_type] = species_adjusted_score
        phenotype_adjusted_scores[phenotype] = embryos_adjusted_scores
    
    return phenotype_adjusted_scores
            



## divide that value by each other
## save to new dictionary with phenotype:embryo_type:phenotype_score