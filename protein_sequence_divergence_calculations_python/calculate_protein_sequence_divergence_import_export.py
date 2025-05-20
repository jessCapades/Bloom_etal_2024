import glob
import os
import transcript_paring
from Bio import SeqIO
import pandas as pd
import sys 

## import define user input files
def get_file_path():
    if len(sys.argv) < 2:
        print("Please provide a sequences file directory, e.g. /path/to/directory")
        exit()
    else:
        file_info = {"sequence_file_directory": sys.argv[1]}
        file_path = file_info["sequence_file_directory"]

        if os.path.exists(file_path):
            return file_info
        else:
            print(f"Unable to access file at {file_path}")
            exit()

## create new directory for files with selected species sequences
def create_dir_if_not_exists(directory_path):
    """Creates a directory if it doesn't exist."""
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

def read_fasta(fasta_file):
    # Parse the FASTA file using Biopython
    sequences = []
    sequence_names = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))  # Convert to string for easier processing
        sequence_names.append(record.id)
    return sequences, sequence_names

## get sequence file locations
def extract_sequence_files(sequence_directory):
    sequence_files = glob.glob(f'{sequence_directory}/*')

    return sequence_files

## create new fasta file name for files with selected species sequences
def create_selected_sequences_fasta_file_names(sequence_file, selected_sequences_directory):
    
    basic_sequence_filename = os.path.basename(os.path.normpath(sequence_file))
    selected_sequence_fasta_filename = f'{basic_sequence_filename}_selected_species.txt'
    selected_sequences_fasta_file_path = f'{selected_sequences_directory}/{selected_sequence_fasta_filename}'

    return selected_sequences_fasta_file_path
    
## get sequence name associated with each sequence
def parse_sequences_file(sequences_file):
    sequences = {}

    with open(sequences_file) as file:
        file = map(str.rstrip, file)

        for line in file:
            if line.startswith('>'):
                sequence_name = line.lstrip('>')
            else:
                sequence_data = sequences.get(sequence_name, '')
                sequence_data += line
                sequences.update({ sequence_name: sequence_data })

    return sequences

def get_sequence_names(fasta_file):
    """Gets a list of sequence names from a FASTA file."""

    sequence_names = []
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence_names.append(record.id)
    return sequence_names

# selects all trancripts associated with species name from fasta file
def filter_sequences(supergroup_sequence_names, longest_sequences):
    filtered_sequences = {}

    for sequence_name, sequence in longest_sequences.items():
        sequence_prefix = transcript_paring.get_sequence_prefix(sequence_name)

        if sequence_prefix in supergroup_sequence_names:
            if sequence_name in filtered_sequences:
                filtered_sequences[sequence_name].append(sequence)
            else:
                filtered_sequences.update({sequence_name: sequence})


    return filtered_sequences

## get species names from input file
def parse_species_names(species_file_name):
    supergroup_sequence_names = set()

    with open(species_file_name) as file:
        file = map(str.rstrip, file)

        for line in file:
            supergroup_sequence_names.add(line)

    return supergroup_sequence_names

## associate all transcript names with the species that are from
def sort_transcript_names_by_species(parsed_sequences):
    transcripts_from_each_species = {}

    for sequence_name in parsed_sequences:

        sequence_prefix = transcript_paring.get_sequence_prefix(sequence_name)
        if sequence_prefix in transcripts_from_each_species:
            transcripts = transcripts_from_each_species[sequence_prefix]
            if isinstance(transcripts, list):
                transcripts.append(sequence_name)
                transcripts_from_each_species.update({sequence_prefix:transcripts})
            else:
                all_transcripts = []
                all_transcripts.append(transcripts)
                transcripts_from_each_species.update({sequence_prefix:all_transcripts})
        else:
            all_transcripts = []
            all_transcripts.append(sequence_name)
            transcripts_from_each_species.update({sequence_prefix:all_transcripts})
    
    return transcripts_from_each_species

def create_fasta_file(file_location, sequences):
    with open(file_location, 'w') as f:
        for k, v in sequences.items():
            f.write('>%s\n%s\n' % (k, v))

def get_filename(filename):
    base, ext = os.path.splitext(os.path.basename(filename))
    return base

def export_dictionary_to_excel_file (species_name_dictionary, file_path, filename):
    # Convert dictionary to a DataFrame
    species_names_and_div_scores = pd.DataFrame(species_name_dictionary)

    # Export DataFrame to Excel
    species_names_and_div_scores.to_excel(f'{file_path}/{filename}')

## create a fasta file of only sequences of species selected in user input file
def create_fasta_of_selected_species(sequences_file, selected_seq_directory, species_name_file):
    
    sequences = parse_sequences_file(sequences_file)
    selected_species_sequence_names = parse_species_names(species_name_file)
    selected_species_sequences = filter_sequences(selected_species_sequence_names, sequences)
    fasta_of_selected_species = create_selected_sequences_fasta_file_names(sequences_file, selected_seq_directory)

    create_fasta_file(fasta_of_selected_species, selected_species_sequences)

def sort_sequences_by_species_name(sequences_file):

    sequences = parse_sequences_file(sequences_file)
    sorted_species_transcripts = sort_transcript_names_by_species(sequences)

    return sequences, sorted_species_transcripts

def create_new_graph_file_names(sequence_file, selected_sequences_directory):
    
    basic_sequence_filename = os.path.basename(os.path.normpath(sequence_file))
    selected_sequence_fasta_filename = f'{basic_sequence_filename}'
    selected_sequences_graph_filepath = f'{selected_sequences_directory}/{selected_sequence_fasta_filename}'

    return selected_sequences_graph_filepath, selected_sequence_fasta_filename