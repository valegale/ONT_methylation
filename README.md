# Nextflow Pipeline for Methylated Motif Extraction and Statistical Analysis from ONT data

This pipeline processes basecalled reads generated by **Dorado**, with basecalled DNA modifications, to identify specific motifs using **Modkit**. Additionally, it performs basic statistical analyses regarding the methylation status of the genome and the extracted motifs.

## Inputs

- **Basecalled reads**: Output from Dorado, including basecalled DNA modifications.
- **Reference file**: A reference genome or sequence against which the reads will be aligned.

## Outputs

- **Motif extraction**: Identified motifs based on DNA modifications using Modkit.
- **Statistical analysis**: Summary statistics of the methylation status of the sample and the extracted motifs.

## Requirements

- **Nextflow**

## How to Run

1. Clone this repository.
2. Basecall your reads with Dorado.
   ```bash
   dorado basecaller sup,6mA,4mC_5mC <pod5_folder> > results.bam

3. Run the pipeline using Nextflow