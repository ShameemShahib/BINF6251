# Implementation Progress Report

## Project Snapshot

**Research Question**
If models are provided with clear descriptions of regulation and protein function, do predictions improve beyond insights derived exclusively from RNA?

**Algorithm**
This project uses a standard Hidden Markov Model (HMM) to summarize DNA methylation data into per-patient regulatory features, a profile HMM to score somatic mutations against conserved protein domains, and a set of classifiers (logistic regression, random forest, neural network) to predict PAM50 breast cancer subtype from RNA expression alone versus RNA combined with those summaries.

**Current Status**
Data loading and the standard HMM are fully implemented and validated on a prototype cohort. The profile HMM and classification stages are yet to be implemented

---

## What Is Implemented

Data Loading - Data_Loading.py
The script loads all four data types directly from UCSC Xena Browser URLs: RNA-seq, DNA methylation (450k), somatic mutations (MC3), and clinical annotations. It intersects patients across RNA, methylation, and PAM50 labels (mutations are not required for intersection since zero mutations is valid data), producing a final cohort of 620 patients. The processed files are saved locally to a Data directory.

HMM - HMM.py
The script implements a HMM from using Baum-Welch for training and Viterbi for decoding. It filters the methylation matrix to the 1,000 most variable CpG sites for the prototype data, trains on a subsample of patients and outputs a single per-patient fraction of active sites dataset that is saved to the Data directory

Not Yet Implemented
- Profile HMM for somatic mutation domain scoring
- Feature matrix construction across all four model configurations
- Classifier training and cross-validated evaluation

---

## Prototype Demo Description

Scripts to run, in order:
1. `Data_Loading.py` — loads and aligns all four data types, saves to Data directory
2. `HMM.py` — trains HMM and decodes patients, saves regulatory features to Data directory

Input files expected:
All data is loaded directly from public URLs defined at the top of Data_Loading.py. No manual downloads are required. After the first run, processed files are saved locally to the Data directory and can be reloaded from there.

Output
HMM.py prints a summary of the fraction of active sites feature across all prototype patients and saves the full results to Data/regulatory_features.tsv

**Example output from prototype run (50 patients, 1,000 CpG sites):**
```
count    50.000000
mean      0.447037
std       0.168354
min       0.126126
25%       0.358358
50%       0.421131
75%       0.577010
max       0.748497
Name: frac_active, dtype: float64
```

---

## Data Documentation

**Origin**
All data comes from The Cancer Genome Atlas Breast Cancer (TCGA-BRCA) cohort, accessed through the UCSC Xena Browser. The four datasets used are:
 - gene expression RNAseq (IlluminaHiSeq, n=1,218)
 - DNA methylation 450k (n=888)
 - somatic mutations MC3 public version (n=791)
 - clinical phenotypes including PAM50 subtype labels (n=1,247).

**Preprocessing**
Patient barcodes are in a consistent 15-character TCGA format across all four files, so no truncation was needed. Patients are intersected across RNA, methylation, and PAM50 labels, yielding 620 patients. Mutations are not required for intersection, patients with no mutations receive a zero feature vector at the domain scoring stage. For the HMM, the methylation matrix is filtered to the 1,000 most variable CpG sites in prototype mode and the first 50 patients are used to keep runtime short.

**Ground Truth**
PAM50 subtype labels from the clinical annotations file serve as ground truth for classification. These labels (LumA, LumB, Basal, HER2, Normal) were derived from RNA-seq data by the original TCGA study and are the standard reference used across the breast cancer literature.

---

## Initial Observations

The prototype HMM output looks biologically reasonable. The mean `frac_active` of 0.45 indicates the model is splitting sites roughly evenly between active and repressed states. The spread from 0.13 to 0.74 shows patient-level variation in regulatory activity.If the HMM were broken, all patients would return identical or skewed values.

One runtime observation was that the forward and backward passes are implemented as loops over CpG sites, which is readable but slow on large sequences. On the prototype with only 1,000 sites and 20 training patients the runtime is only a few minutes. However, scaling to 10,000 sites and 200 patients for the full run will be significantly slower and may require optimization.

---

## Reflection on Changes and Challenges

**Deviation from pseudocode**
The original pseudocode specified three per-patient regulatory features: frac_active, frac_repressed, and n_transitions. The current implementation outputs only frac_active. Since frac_repressed is always 1 - frac_active, it was removed to keep the code concise. n_transitions may be added back later if additional regulatory signal is needed, but was dropped for now to keep the implementation simple.

**Patient cohort size**
The original proposal aimed to tackle around 1000 patients. The intersection of RNA, methylation, and PAM50 labels yielded far less, primarily because not all patients have all the necessary types of data publically available. This is still a workable cohort size and does not affect the project.

**Data access**
Rather than downloading files manually, data is loaded directly from UCSC Xena URLs, which simplified the data acquisition step significantly.

**Barcode matching**
The original plan anticipated needing to truncate TCGA barcodes to 12 characters for matching across files. Thankfully the barcodes were consistent across all the files, so no truncation was needed.

---

## Next Steps

- Implement the profile HMM for somatic mutation domain scoring, producing a max_disruption feature per patient
- Build the four feature matrix configurations: RNA-only, RNA + regulatory, RNA + domain, and full
- Implement stratified 5-fold cross-validation for logistic regression, random forest, and a small neural network
- Compare macro F1 scores across configurations to evaluate whether regulatory and domain features add value beyond RNA alone
- Validate the HMM by checking that patients with higher frac_active scores tend to have higher overall gene expression
