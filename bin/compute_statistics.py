#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
import argparse
import os

def read_modkit(modkit_output):
    # Valid coverage is = Nmod + Nother_mod + Ncanonical.  In modkit the value: fraction_modified is the one in the bedgraph, but here
    # we compute percent_modified as a fraction of the modified bases / total_coverage, rather than valid coverage. 
    # Total coverage include also the modification below threshold and other bases, so it is more accurate. 

    d = pd.read_csv(modkit_output, sep="\t", header=None)

    d.columns = [
        "Contig", "Position", "End", "Modification", "drop0", "Strand", "drop1", "drop2", "drop3", "Valid_coverage",  
        "fraction_modified", "Modified_bases", "Unmodified_bases", "Other_mod_base", "drop6", "Modification_below_threshold", "Other_bases", "drop7"]
    
    d = d[d.columns[~d.columns.str.contains("drop")]]
    
    d["Total_coverage"] = (
        d["Modified_bases"]
        + d["Unmodified_bases"]
        + d["Other_mod_base"]
        + d["Modification_below_threshold"]
        + d["Other_bases"]
    )

    d["Percent_modified"] = d["Modified_bases"] / d["Total_coverage"]  
    d["fraction_modified"] = d["fraction_modified"] / 100  # this is the value computed by modkit (-bedgraph)
    
    return d

def compute_statistics_modifications(filtered_table, results_folder):
    # % of methylated a (6mA)
    a_modifications = filtered_table[filtered_table["Modification"] == 'a']
        
    a_modification_counts = a_modifications["Contig"].value_counts().reset_index()
    a_modification_counts.columns = ['Contig', 'tot_6mA']
    a_modification_counts['tot_a'] = a_modification_counts['Contig'].map(lambda x: reference_genome[x].lower().count('a') if x in reference_genome else None)
    a_modification_counts['%_methylation'] = (a_modification_counts['tot_6mA'] / a_modification_counts['tot_a']) * 100
        
    # % of methylated c (5mC)
    c5_modifications = filtered_table[filtered_table["Modification"] == 'm']
    c5_modification_counts = c5_modifications["Contig"].value_counts().reset_index()
    c5_modification_counts.columns = ['Contig', 'tot_5mC']
    c5_modification_counts['tot_c'] = c5_modification_counts['Contig'].map(lambda x: reference_genome[x].lower().count('c') if x in reference_genome else None)
    c5_modification_counts['%_methylation'] = (c5_modification_counts['tot_5mC'] / c5_modification_counts['tot_c']) * 100

        
    # % of methylated c (4mC)
    c4_modifications = filtered_table[filtered_table["Modification"] == '21839']
    c4_modification_counts = c4_modifications["Contig"].value_counts().reset_index()
    c4_modification_counts.columns = ['Contig', 'tot_4mC']
    c4_modification_counts['tot_c'] = c4_modification_counts['Contig'].map(lambda x: reference_genome[x].lower().count('c') if x in reference_genome else None)
    c4_modification_counts['%_methylation'] = (c4_modification_counts['tot_4mC'] / c4_modification_counts['tot_c']) * 100
 

    output_path =  os.path.join(results_folder, "6mA_percentage.csv")
    a_modification_counts.to_csv(output_path, sep="\t", index = False)
    output_path =  os.path.join(results_folder, "5mC_percentage.csv")
    c5_modification_counts.to_csv(output_path, sep="\t", index = False)
    output_path =  os.path.join(results_folder, "4mC_percentage.csv")
    c4_modification_counts.to_csv(output_path, sep="\t", index = False)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Read Modkit output and process the table')
    parser.add_argument('modkit_output_path', type=str, help='Path to the Modkit output file')
    parser.add_argument('fasta_file_path', type=str, help='Path to the reference genome FASTA file')
    parser.add_argument('results_folder', type=str, help='Name of the folder with the results')
    parser.add_argument('--min_coverage', type=int, default=10, help='The minimum coverage required to consider a site for methylation (default 10)')
    parser.add_argument('--percent_cutoff', type=int, default=0.5, help='Filter threshold for positions with a low level of methylation (default 0.5)') 
    args = parser.parse_args()


    modkit_file_path = args.modkit_output_path
    fasta_file_path = args.fasta_file_path
    min_coverage = args.min_coverage
    percent_cutoff = args.percent_cutoff

    if not os.path.exists(args.results_folder):
        os.makedirs(args.results_folder)        

    reference_genome = {}
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        reference_genome[record.id] = record.seq


    modkit_table = read_modkit(modkit_file_path)
    filtered_table = modkit_table[(modkit_table.Total_coverage >= min_coverage) & (modkit_table.Percent_modified >= percent_cutoff)].drop("End", axis=1) 
    compute_statistics_modifications(filtered_table, args.results_folder)