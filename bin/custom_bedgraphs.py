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
        "Fraction_modified", "Modified_bases", "Unmodified_bases", "Other_mod_base", "drop6", "Modification_below_threshold", "Other_bases", "drop7"]
    
    d = d[d.columns[~d.columns.str.contains("drop")]]
    
    d["Total_coverage"] = (
        d["Modified_bases"]
        + d["Unmodified_bases"]
        + d["Other_mod_base"]
        + d["Modification_below_threshold"]
        + d["Other_bases"]
    )

    d["Percent_modified"] = d["Modified_bases"] / d["Total_coverage"]  
    d["Fraction_modified"] = d["Fraction_modified"] / 100  # this is the value computed by modkit (-bedgraph)
    
    return d


def save_filtered_tables(modkit_table, folder_tables, percent_cutoff, min_coverage):
    modifications = ["a", "m", "21839"] 
    modification_names = {"a": "6mA", "m": "5mC", "21839": "4mC"}

    for modification in modifications: 
        filtered_table = modkit_table[(modkit_table.Total_coverage >= min_coverage) & (modkit_table.Modification == modification) & (modkit_table.Percent_modified >= percent_cutoff)].drop("End", axis=1)   
        filtered_table_file_path  = f"{folder_tables}/filtered_modkit_table_{modification_names[modification]}.tsv"
        filtered_table.to_csv(filtered_table_file_path, sep='\t', index=False)


def create_bedgraphs_file(modkit_table, folder_bedgraphs, prefix, min_coverage):
    strands = ["-", "+"]
    modifications = ["a", "m", "21839"]
    columns_bedgraph = ['Contig', 'Position', 'End', 'Percent_modified', 'Total_coverage']

    for modification in modifications:
        for strand in strands:
            bedgraph_table = modkit_table[(modkit_table.Strand == strand) & (modkit_table.Modification == modification) & (modkit_table.Total_coverage >= min_coverage)][columns_bedgraph]

            file_name = f"{prefix}_{'positive_' if strand == '+' else 'negative_'}{modification}.bedgraph"
            file_path = f"{folder_bedgraphs}/{file_name}"

            bedgraph_table.to_csv(file_path, sep='\t', header=False, index=False)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Read Modkit output and process the table')
    parser.add_argument('modkit_output_path', type=str, help='Path to the Modkit output file')
    parser.add_argument('fasta_file_path', type=str, help='Path to the reference genome FASTA file')
    parser.add_argument('results_folder', type=str, help='Name of the folder with the results')
    parser.add_argument('--min_coverage', type=int, default=10, help='The minimum coverage required to consider a site for methylation (default 10)')
    parser.add_argument('--percent_cutoff', type=float, default=0, help='Filter threshold for positions with a low level of methylation (default 0)') 
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
    

    results_table_path =  os.path.join(args.results_folder, "modifications_tables")
    save_filtered_tables(modkit_table, results_table_path, percent_cutoff, min_coverage)

    results_bedgraphs_path =  os.path.join(args.results_folder, "bedgraphs_customized")
    create_bedgraphs_file(modkit_table, results_bedgraphs_path, "bedgraph", min_coverage)