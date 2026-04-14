# BINF6251
Project for BINF6251 - Lab for Algorithms

# Multi-Omics Breast Cancer Subtype Classification

## Project Overview

### Research Question
If models are provided with clear descriptions of regulation and protein 
function, do predictions improve beyond insights derived exclusively from RNA?

### Algorithm and Approach
This project implements three algorithms from scratch or near-scratch:

1. **Standard Hidden Markov Model**
   Processes DNA methylation data to infer regulatory state (active vs 
   repressed) at each genomic position per patient, producing three 
   summary features: fraction of active sites, fraction of repressed 
   sites, and number of regulatory state transitions.

2. **Profile HMM via PyHMMER + PolyPhen scores**  
   Processes somatic mutation data to score whether mutations fall inside 
   known Pfam protein domains, producing three disruption features per 
   patient: maximum disruption score, number of domain hits, and mean 
   disruption score.

3. **Neural Network Classifier** 
   Trains a two-layer feedforward network to predict PAM50 breast cancer 
   subtype (LumA, LumB, Her2, Basal, Normal) across four feature 
   configurations: RNA-only, RNA + regulatory, RNA + domain, and full 
   multi-omics.

### Data
- **Input:** TCGA-BRCA RNA expression, DNA methylation (450k array), 
  somatic mutations (MAF format), and clinical annotations (PAM50 labels)
- **Source:** UCSC Xena Browser (TCGA-BRCA cohort)
- **Output:** Per-patient macro F1 scores across four feature 
  configurations, and per-patient regulatory and domain feature vectors

---

## Installation and Setup

### Requirements
- Python 3.10 or higher

Step 1: Clone the repository
Step 2: Create and activate a virtual environment
Step 3: Install dependencies in requirements.txt
Step 4: Download Pfam database (required for profile HMM)
```bash
curl -O https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
```
Step 5: Run Data_Loading.py to create a data directory with all the necessary data (over 2GB)

## Quick Start

After completing installation and data setup, run the three scripts in order:
HMM.py -> ProfileHMM.py -> NeuralNetwork.py

Step 1: Train HMM and generate regulatory features
- The output should be a file named regulatory_features.tsv in the Data directory

Step 2: Generate domain disruption features using the profile HMM
- The output should be a file named domain_features.tsv in the Data directory

Step 3: Train the classifier and review results
- The output will be printed and include a results table

---

## Usage and Options

### HMM.py

Trains an HMM on DNA methylation data using Baum-Welch, then decodes all patients with Viterbi. The trained model is saved as a .pkl file and can be reloaded without retraining:

### ProfileHMM.py

Fetches protein sequences, scans against Pfam domains, and scores somatic mutations by domain membership and PolyPhen pathogenicity, which predicts how harmful a mutation in an amino acid can be for a protein's overall structure and function. Protein and domain caches are saved as .pkl files as well.

### NeuralNetwork.py

Trains a two-layer neural network across four feature configurations using 5-fold stratified cross-validation.

## Limitations and Assumptions

One of the biggest limitations with handling data of this scale is computational power. The pipeline was developed and tested on a MacBook Air M4 with 16GB of RAM, and the full workflow took over an hour to complete. All results are based on a single training run using a fixed random seed (42), so there may be slight variation if the model is rerun with a different seed.

On the methylation side, the HMM assumes that consecutive CpG sites are genomic neighbors. Because of this, sites need to be sorted by chromosomal position before training, and any sites on non-standard contigs are excluded to keep the data consistent.

On the mutation side, many patients have no domain-disrupting mutations, which often results in sparse or entirely zero feature vectors and reduces the ability of these features to effectively distinguish between samples.

In addition, the MAF file used for mutation data may include mutations from multiple tumor types. To make the analysis more relevant, strict filtering was applied, and only genes mutated in at least two patients were kept. 
