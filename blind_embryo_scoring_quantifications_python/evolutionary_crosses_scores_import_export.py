import pandas as pd
import sys
import os

def get_file_path():
    if len(sys.argv) < 3:
        print("Please provide a directory and file path, e.g. /path/to/directory filename.ext")
        exit()
    else:
        file_info = {"directory": sys.argv[1], "filename": sys.argv[2]}
        file_path = file_info["directory"] + file_info["filename"]

        if os.path.exists(file_path):
            return file_info
        else:
            print(f"Unable to access file at {file_path}")
            exit()

def import_phenotype_scores(phenotype_scores_directory, phenotype_scores_filename, embryos_for_exclusion):
    
    file_path = f"{phenotype_scores_directory}/{phenotype_scores_filename}"
    
    phenotype_scores = pd.read_csv(file_path)
    phenotype_scores = phenotype_scores[~(phenotype_scores["embryo_type"].isin(embryos_for_exclusion))]
 
    
    return phenotype_scores

def export_adjusted_phenotype_scores(adjusted_phenotype_scores, file_path):

    phenotype_scores_for_export = pd.DataFrame.from_dict(adjusted_phenotype_scores, orient='index')

    phenotype_scores_for_export.to_csv(file_path)


def import_adjusted_phenotype_scores(phenotype_scores_file_path):

    adjusted_phenotype_scores = pd.read_csv(phenotype_scores_file_path)
    
    return adjusted_phenotype_scores