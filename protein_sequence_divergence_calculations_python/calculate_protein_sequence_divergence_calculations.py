from Bio import SeqIO
from Bio.Align import substitution_matrices
from Bio import Align
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
import calculate_protein_sequence_divergence_import_export as imp_exp
import os
import transcript_paring as trans_paring
# from Bio import Phylo
# from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
# from Bio.Phylo.TreeConstruction import _DistanceMatrix
# from ete3 import Tree

def define_cel_ortholog(sequence_file_directory):
    cel_transcript_names = []
    sequence_file_names = os.listdir(sequence_file_directory)
    
    for filename in sequence_file_names:
        cel_transcript_name = filename.split("_")[0]
        cel_transcript_names.append(cel_transcript_name)

    return cel_transcript_names

def retreive_matching_cel_transcript_to_current_file(cel_transcript_names, sequence_file_path):

    for transcript_name in cel_transcript_names:
        if transcript_name in sequence_file_path:
            cel_transcript_name = f'CELEG.{transcript_name}'
            print(cel_transcript_name)
        else:
            print("not this transcript")
    
    return cel_transcript_name

def retrieve_cel_sequence (full_cel_transcript_name, transcript_name_to_sequences):

    if any(key.startswith(full_cel_transcript_name) for key in transcript_name_to_sequences):
        cel_sequence_for_comparison = next(v for k,v in transcript_name_to_sequences.items() if full_cel_transcript_name in k)
    else:
        print(f'{full_cel_transcript_name} thats not a real gene')
        cel_sequence_for_comparison = []
    
    return cel_sequence_for_comparison

def identify_and_pull_cel_ortholog(cel_short_transcript_names, transcript_names_and_sequences, sequence_file_path):
    cel_ortholog_name = retreive_matching_cel_transcript_to_current_file(cel_short_transcript_names, sequence_file_path)
    cel_ortholog_for_pruning = retrieve_cel_sequence(cel_ortholog_name, transcript_names_and_sequences)

    return cel_ortholog_name, cel_ortholog_for_pruning

# calculate divergence scores only comparing to C. elegans transcript
def calculate_blosum_divergence_for_cel(sequences, cel_ortholog_for_pruning):
    blosum62 = substitution_matrices.load("BLOSUM62")
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.substitution_matrix = blosum62
    aligner.open_gap_score = -12
    aligner.extend_gap_score = -0.5
    number_of_sequences = len(sequences)
    distance_matrix = np.zeros(number_of_sequences)

    for i in range(number_of_sequences):
        # Perform global sequence alignment using BLOSUM62
        alignments = aligner.align(sequences[i], cel_ortholog_for_pruning)
        if len(alignments) > 10000:
            number_of_alignments = 10000
        else:
            number_of_alignments = len(alignments)
        best_alignment_score = np.max(np.array([alignments[k].score for k in range(number_of_alignments)]))

        distance_matrix[i] = best_alignment_score
    
    max_possible_score = np.min(np.array(distance_matrix)) # Approximate max score for BLOSUM62
    for i in range(number_of_sequences):
        # Normalize score by the maximize score for each sequence pair
        similarity = distance_matrix[i] / max_possible_score
        divergence = 1 - similarity

        # Store the divergence (distance)
        distance_matrix[i] = divergence

    return distance_matrix

def collect_sequence_divergence_scores(cel_distance_matrix, transcript_names):
    transcript_names_and_divergence_scores = {}

    if len(cel_distance_matrix) == len(transcript_names):
        
        for i in range(len(transcript_names)):
            transcript_name = transcript_names[i]
            divergence_score = cel_distance_matrix[i]
            transcript_names_and_divergence_scores.update({transcript_name: divergence_score})
    else:
        print("sequence name list length doesn't match divergence matrix length")
    
    return transcript_names_and_divergence_scores

def subset_dictionary_based_on_species(transcript_according_to_species, transcript_divergence_scores):
    
    best_orthologs = []
    
    for species_name in transcript_according_to_species:
        transcript_names = transcript_according_to_species[species_name]
        species_specific_divergence_scores = {k: transcript_divergence_scores[k] for k in (transcript_names)}

        if len(species_specific_divergence_scores) == 1: 

            for transcript_name in species_specific_divergence_scores:
                best_orthologs.append(transcript_name)
        else:
            most_similar_ortholog = min(species_specific_divergence_scores, key=species_specific_divergence_scores.get)
            best_orthologs.append(most_similar_ortholog)
    
    return best_orthologs

def select_sequences_based_on_highest_alignment_score(best_orthologs, transcripts_and_sequences):

    sequences_from_best_orthologs = {k: transcripts_and_sequences[k] for k in (best_orthologs)}

    return sequences_from_best_orthologs



def select_best_orthologs_based_on_distance_scores (celegans_distance_matrix, transcript_names, transcripts_under_species_names, transcript_names_and_sequences):

    transcripts_and_divergence_scores = collect_sequence_divergence_scores(celegans_distance_matrix, transcript_names)
    best_orthologs_per_species = subset_dictionary_based_on_species(transcripts_under_species_names, transcripts_and_divergence_scores)
    sequences_from_best_orthologs = select_sequences_based_on_highest_alignment_score(best_orthologs_per_species, transcript_names_and_sequences)

    return sequences_from_best_orthologs

# Function to calculate divergence using BLOSUM62 scoring
def calculate_blosum_divergence(sequences):
    blosum62 = substitution_matrices.load("BLOSUM62")
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.substitution_matrix = blosum62
    number_of_sequences = len(sequences)
    aligner.open_gap_score = -12
    aligner.extend_gap_score = -0.5
    distance_matrix = np.zeros((number_of_sequences, number_of_sequences))

    for i in range(number_of_sequences):
        for j in range(i, number_of_sequences):
            # Perform global sequence alignment using BLOSUM62
            alignments = aligner.align(sequences[i], sequences[j])
            if len(alignments) > 10000:
                number_of_alignments = 10000
            else:
                number_of_alignments = len(alignments)
            best_alignment_score = np.max(np.array([alignments[k].score for k in range(number_of_alignments)]))
            distance_matrix[i][j] = best_alignment_score
    
    for i in range(number_of_sequences):
        for j in range(i+1, number_of_sequences):
            # Normalize score by the maximize score for each sequence pair
            max_possible_score = np.min(np.array([distance_matrix[i][i], distance_matrix[j][j]]))  # Approximate max score for BLOSUM62
            similarity = distance_matrix[i][j] / max_possible_score
            divergence = 1 - similarity

            # Store the divergence (distance)
            distance_matrix[i][j] = divergence
            distance_matrix[j][i] = divergence
   # Replace self self matches with 0
    for i in range(number_of_sequences):
        distance_matrix[i][i] = 0

    return distance_matrix

# Step 3: Visualize divergence using hierarchical clustering
def plot_dendrogram(distance_matrix, labels, sequence_filename):
    # Perform hierarchical clustering
    converted_1D_distance_matrix = squareform(distance_matrix)
    linked = linkage(converted_1D_distance_matrix, 'average')

    # Create a dendrogram
    plt.figure(figsize=(8, 6))
    dendrogram(linked, labels=labels, orientation='top')
    plt.ylim(0, 1.2)
    plt.title(f'{sequence_filename}')
    plt.xlabel('Sequence')
    plt.ylabel('Divergence')
    

def create_divergence_heatmap(distance_matrix, species_name_labels, filename):

    # Create the heatmap using imshow
    plt.figure(figsize=(8, 6))  # Optional: set figure size
    plt.imshow(distance_matrix, cmap='coolwarm', interpolation='nearest', vmin=0, vmax=1.0)
    # Add color bar to the side
    plt.colorbar()
    # Optional: Add labels and title
    plt.title(f'{filename}', fontsize=16)
    plt.xlabel("Species", fontsize=12)
    plt.ylabel("Species", fontsize=12)
    # Optional: Set ticks for rows and columns
    plt.xticks(ticks=np.arange(distance_matrix.shape[1]), labels=species_name_labels)
    plt.yticks(ticks=np.arange(distance_matrix.shape[0]), labels=species_name_labels)


def standardize_species_names(pruned_sequence_filename_and_transcripts):

    for gene_name in pruned_sequence_filename_and_transcripts:

        transcript_names = pruned_sequence_filename_and_transcripts[gene_name]
        transcript_prefixes = []

        for transcript in transcript_names:
            species_name_prefix = trans_paring.get_sequence_prefix(transcript)
            transcript_prefixes.append(species_name_prefix)
        pruned_sequence_filename_and_transcripts[gene_name] = transcript_prefixes
    
    return pruned_sequence_filename_and_transcripts

def collect_cbn_divergence_scores(genes_names_and_standard_species_names, genes_names_and_divergence_scores):

    gene_names_and_cbrenneri_scores = {}
    
    for gene_name in genes_names_and_standard_species_names:
        species_names = genes_names_and_standard_species_names[gene_name]
        # associate each score with the species name
        # associate each species name with the index from the divergence matrix
        if "CBREN" in species_names:
            distance_matrix = genes_names_and_divergence_scores[gene_name]
            cbrenneri_divergence_index = species_names.index("CBREN") 
            cbrenneri_divergence_scores = distance_matrix[cbrenneri_divergence_index]
            gene_names_and_cbrenneri_scores[gene_name] = cbrenneri_divergence_scores
        else:
            print("CBREN does not have an ortholog")
    
    return gene_names_and_cbrenneri_scores
        

def remove_genes_without_orthologs_for_all_species(cbn_divergence_scores):

    cbn_scores_versus_all_species = {}

    for gene_name in cbn_divergence_scores:
        divergence_scores = cbn_divergence_scores[gene_name]
        if len(divergence_scores) < 5:
            print("not all orthologs present for this gene")
        else:
            cbn_scores_versus_all_species[gene_name] = divergence_scores
    
    return cbn_scores_versus_all_species

def associate_species_names_with_transcript_scores(gene_names_and_standard_species_names, gene_names_with_divergence_scores):
    
    cbn_vs_species_associated_scores = {}

    for gene_name in gene_names_and_standard_species_names:
        species_and_divergence_score = {}

        species_names = gene_names_and_standard_species_names[gene_name]
        for i in range(len(species_names)):

            divergence_score = gene_names_with_divergence_scores[gene_name][i]
            species_and_divergence_score [species_names[i]] = divergence_score
        cbn_vs_species_associated_scores [gene_name] = species_and_divergence_score

    return cbn_vs_species_associated_scores 

