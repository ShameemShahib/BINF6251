import os
import pandas as pd

# obtain URL for each data file from UCSC Xena Browser
rna_url = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2.gz"
methyl_url = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHumanMethylation450.gz"
mutation_url = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/mc3%2FBRCA_mc3.txt.gz"
clin_annotation_url = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix"

# Load the rna, methylation, mutation, and clinical annotation data
rna = pd.read_csv(rna_url, sep='\t', index_col=0, compression='gzip')
methyl = pd.read_csv(methyl_url, sep='\t', index_col=0, compression='gzip')
mutation = pd.read_csv(mutation_url, sep='\t', compression='gzip')
clin_annotation = pd.read_csv(clin_annotation_url, sep='\t', index_col=0)

# find patients with PAM50 labels(cancer subtypes) and ignore the rest
PAM_50 = clin_annotation['PAM50Call_RNAseq'].dropna()

# find patients that intersect with rna, methylation, and PAM50 data
Intersect_cohort = sorted(set(rna.columns) & set(methyl.columns) & set(PAM_50.index))

# subset RNA, methyl, PAM50
# for mutations, even 0 mutations is useful data, so include all patients who fit the intersected data cohort
rna = rna[Intersect_cohort]
methyl = methyl[Intersect_cohort]
mutation = mutation[mutation['sample'].isin(Intersect_cohort)]
PAM_50 = PAM_50[Intersect_cohort]

# Save reduced dataset into new directory
os.makedirs('Data', exist_ok=True)
rna.to_csv('Data/RNAseq.tsv', sep='\t')
methyl.to_csv('Data/methyl.tsv', sep='\t')
mutation.to_csv('Data/mutation.tsv', sep='\t', index=False)
clin_annotation.to_csv('Data/clinicalMatrix.tsv', sep='\t')
