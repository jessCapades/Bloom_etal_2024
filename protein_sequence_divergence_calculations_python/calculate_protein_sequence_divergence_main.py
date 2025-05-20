import calculate_protein_sequence_divergence_import_export as imp_exp
import calculate_protein_sequence_divergence_calculations as calculations
import os
import matplotlib.pyplot as plt
import numpy as np

# collect file info and file paths
sequence_evolution_file_info = imp_exp.get_file_path()
sequence_evolution_directory = sequence_evolution_file_info["sequence_file_directory"]
sequence_files_directory = f'{sequence_evolution_directory}FASTANexusForJess/Prots'
parsed_sequence_files_directory = f'{sequence_evolution_directory}bloom_etal_sequences_for_divergence_analysis'
selected_species_filename = f'{sequence_evolution_directory}bloom_etal_2025_hyb_species.txt'
divergence_analysis_results = f'{sequence_evolution_directory}bloom_etal_divergence_analysis_results'
divergence_analysis_results_dendrograms = f'{divergence_analysis_results}/dendrograms'
divergence_analysis_results_heatmaps = f'{divergence_analysis_results}/heatmaps'
pruned_sequences_for_final_alignment = f'{sequence_evolution_directory}/bloom_etal_divergence_analysis_pruned_orthologs'

imp_exp.create_dir_if_not_exists(parsed_sequence_files_directory)
sequences_for_seq_divergence = imp_exp.extract_sequence_files(sequence_files_directory)

for sequence_file in sequences_for_seq_divergence:
    
    parsed_sequence_file = imp_exp.create_fasta_of_selected_species(sequence_file, parsed_sequence_files_directory, selected_species_filename)


selected_species_seq_for_divergence = imp_exp.extract_sequence_files(f'{parsed_sequence_files_directory}')
cel_transcripts = calculations.define_cel_ortholog(parsed_sequence_files_directory)

imp_exp.create_dir_if_not_exists(divergence_analysis_results_dendrograms)
imp_exp.create_dir_if_not_exists(divergence_analysis_results_heatmaps)
imp_exp.create_dir_if_not_exists(pruned_sequences_for_final_alignment)

# determine best orthologs from orthofinder for sequence divergence analysis
for sequence_file in selected_species_seq_for_divergence:

   fasta_file, sequence_names = imp_exp.read_fasta(sequence_file)
   sequence_names = imp_exp.get_sequence_names(sequence_file)
   sequences, transcripts_by_species = imp_exp.sort_sequences_by_species_name(sequence_file)
   cel_transcript_name, cel_ortholog_for_pruning = calculations.identify_and_pull_cel_ortholog(cel_transcripts, sequences, sequence_file)
   cel_divergence_matrix = calculations.calculate_blosum_divergence_for_cel(fasta_file, cel_ortholog_for_pruning)
   best_orthologs_per_species = calculations.select_best_orthologs_based_on_distance_scores(cel_divergence_matrix, sequence_names, transcripts_by_species, sequences)
   
   imp_exp.create_fasta_file(f'{pruned_sequences_for_final_alignment}/{cel_transcript_name}.txt', best_orthologs_per_species)

# pairwise alignment between sequences that align best to the cel ortholog
final_sequences_for_divergence_analysis = imp_exp.extract_sequence_files(f'{pruned_sequences_for_final_alignment}')

pruned_sequence_filenames_and_divergence_scores = {}
pruned_sequence_filenames_and_pruned_sequence_names = {}

for pruned_sequence_file in final_sequences_for_divergence_analysis:
   
   pruned_sequence_filename = imp_exp.get_filename(pruned_sequence_file)
   pruned_fasta_file, pruned_sequence_names = imp_exp.read_fasta(pruned_sequence_file)
   divergence_matrix = calculations.calculate_blosum_divergence(pruned_fasta_file)
   distance_matrix_dendrograms, divergence_matrix_filename = imp_exp.create_new_graph_file_names(pruned_sequence_file, divergence_analysis_results_dendrograms)
   calculations.plot_dendrogram(divergence_matrix, pruned_sequence_names, divergence_matrix_filename)
   plt.savefig(f'{distance_matrix_dendrograms}_dendrogram.pdf', dpi=300, transparent=True)
   plt.close()
   distance_matrix_heatmaps, distance_heatmap_filename = imp_exp.create_new_graph_file_names(pruned_sequence_file, divergence_analysis_results_heatmaps)
   calculations.create_divergence_heatmap(divergence_matrix, pruned_sequence_names, divergence_matrix_filename)
   plt.savefig(f'{distance_matrix_heatmaps}_heatmaps.svg', dpi=300, transparent=True)
   plt.close()
   pruned_sequence_filenames_and_pruned_sequence_names[pruned_sequence_filename] = pruned_sequence_names
   pruned_sequence_filenames_and_divergence_scores[pruned_sequence_filename] = divergence_matrix


pruned_sequence_filenames_and_standard_sequence_names = calculations.standardize_species_names(pruned_sequence_filenames_and_pruned_sequence_names) # standardize species names 
cbrenneri_divergence_scores = calculations.collect_cbn_divergence_scores(pruned_sequence_filenames_and_standard_sequence_names, pruned_sequence_filenames_and_divergence_scores) # collect cbrenneri specific div scores
cbrenneri_with_all_orthologs = calculations.remove_genes_without_orthologs_for_all_species(cbrenneri_divergence_scores)
cbrenneri_species_names_with_all_orthologs = calculations.remove_genes_without_orthologs_for_all_species(pruned_sequence_filenames_and_standard_sequence_names)
species_associated_cbn_divergence_scores = calculations.associate_species_names_with_transcript_scores(cbrenneri_species_names_with_all_orthologs, cbrenneri_with_all_orthologs)
imp_exp.export_dictionary_to_excel_file(species_associated_cbn_divergence_scores, divergence_analysis_results, 'species_associated_divergence_scores.xlsx')

gene_labels = list(cbrenneri_with_all_orthologs.keys())
cbrenneri_with_all_orthologs_converted = np.array(list(cbrenneri_with_all_orthologs.values()))
plt.imshow(cbrenneri_with_all_orthologs_converted, cmap='coolwarm', interpolation='nearest')
plt.xlabel("Species", fontsize=12)
plt.ylabel("Genes", fontsize=12)
plt.xticks(ticks=np.arange(cbrenneri_with_all_orthologs_converted.shape[1]), labels= ["Cbn", "Cel", "Cre", "Csi", "Csp48"])
plt.yticks(ticks=np.arange(cbrenneri_with_all_orthologs_converted.shape[0]), labels=gene_labels)
plt.colorbar()
plt.savefig(f'{divergence_analysis_results_heatmaps}/Cbrenneri_all_genes_heatmap.svg', dpi=300, transparent=True)
plt.close()
# visualize this new heat map