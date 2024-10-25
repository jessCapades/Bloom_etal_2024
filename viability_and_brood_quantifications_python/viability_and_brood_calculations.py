## Import libraries
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys


# 1. Import data
# 2. Analyze data
# 	1. Sum data by day
# 	2. Calculate viability by day
#   3. Calculate brood size by day
# 	4. Plot brood size by day
# 	5. Plot viability by day
# 3. Save plots pdf 


## Import dataset
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
viability_file_info = get_file_path()

viability_and_brood_directory = viability_file_info["directory"]
viability_and_brood_filename = viability_file_info["filename"]

viability_and_brood_dataset = pd.read_csv(viability_and_brood_directory + viability_and_brood_filename)  



## Analyze data
def calculate_embryo_viability(embryo_count, unhatched_count):
	viability = None

	if embryo_count == 0:
		viability = 0.0
	else:
		viability = (embryo_count - unhatched_count) / embryo_count

	return viability


descriptor_columns = ['Species',
                      'plate_number']
embryo_counts_columns = [['day_0_embryo_count',
						  'day_0_unhatched_embryo_count'],
						 ['day_1_embryo_count',
						  'day_1_unhatched_embryo_count'],
						 ['day_2_embryo_count',
						  'day_2_unhatched_embryo_count'],
						 ['day_3_embryo_count',
						  'day_3_unhatched_embryos']]
embryo_counts_columns_flat = [day for daylist in embryo_counts_columns for day in daylist]
species_embryo_count_columns = descriptor_columns + embryo_counts_columns_flat
embryo_counts = viability_and_brood_dataset[species_embryo_count_columns]
embryo_count_by_plate = embryo_counts.groupby(['plate_number']).sum(numeric_only = True)

plate_numbers = embryo_count_by_plate.index
plate_viability_by_day = {}
plate_brood_by_day = {}

for plate_number in plate_numbers:
	plate_data = embryo_count_by_plate.loc[plate_number]
	plate_viability_by_day[plate_number] = []
	plate_brood_by_day[plate_number] = []

	for [count_column, unhatched_column] in embryo_counts_columns:
		embryo_count = plate_data.loc[count_column]
		unhatched_count = plate_data.loc[unhatched_column]
		viability = calculate_embryo_viability(embryo_count, unhatched_count)
		plate_viability_by_day[plate_number].append(viability)
		plate_brood_by_day[plate_number].append(embryo_count)


### Convert data to dataframe
plate_viability_by_day_dataframe = pd.DataFrame.from_dict(plate_viability_by_day)
plate_brood_by_day_dataframe = pd.DataFrame.from_dict(plate_brood_by_day)
viability_brood_quantifications_filename = os.path.splitext(viability_and_brood_filename)[0]


plate_viability_by_day_dataframe.to_csv(f'{viability_and_brood_directory+viability_brood_quantifications_filename}_viability_by_day.csv')
plate_brood_by_day_dataframe.to_csv(f'{viability_and_brood_directory+viability_brood_quantifications_filename}_brood_by_day.csv')
