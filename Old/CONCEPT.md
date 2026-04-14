# Conceptual Progress Report & Pseudocode

## Project Recap
#### Research Question

**If models are provided with clear descriptions of regulation and protein function, do predictions improve beyond insights derived exclusively from RNA?**

#### Summary
This project asks whether adding summaries of gene regulation and protein function improves breast cancer subtype prediction beyond what RNA expression alone can achieve. Using The Cancer Genome Atlas (TCGA) BRCA, dataset DNA methylation data is processed through a standard Hidden Markov Model to infer whether genomic regions are in an active or repressed regulatory state. Additionally, somatic mutation data is processed through a profile Hidden Markov Model to evaluate whether mutations disrupt conserved functional protein domains. These summaries are combined with RNA expression features and used to train classifiers predicting breast cancer subtype, allowing direct comparison between RNA-only and multi-omics models.

## Inputs, Outputs, and Assumptions

#### Inputs
The pipeline takes four patient-matched inputs. The RNA expression matrix is a tab-separated file (genes × patients) containing normalized expression values. The DNA methylation matrix is the same format (CpG sites × patients) containing beta values between 0 and 1, sorted by genomic position before being passed to the HMM. Somatic mutations are provided in MAF format, with one row per mutation describing the gene, protein position, and amino acid change. Clinical annotations provide the PAM50 subtype label for each patient, which has five prediction targets (LumA, LumB, HER2, Basal, Normal) that classifiers are trained to predict.

#### Outputs
The pipeline produces per-patient regulatory feature vectors (frac_active, frac_repressed, n_transitions) and domain mutation vectors (max_disruption, n_domain_hits, mean_disruption), four feature matrices corresponding to the four model configurations (RNA-only, RNA + regulatory, RNA + domain, full), and classification metrics (accuracy, macro F1) for each configuration and classifier type.

#### Assumptions
Patient barcodes are truncated to 12 characters before matching across data types, because different TCGA files use different barcode lengths for the same patient. CpG sites must be in genomic order because the HMM treats consecutive sites as neighbors, so shuffling them would break that relationship. Methylation/beta values are assumed to cluster into two groups: one around low values (active regions) and one around high values (repressed regions), which justifies modeling each state as a bell curve. Two hidden states are assumed to be sufficient,since the model is not trying to capture every possible regulatory nuance, just a broad active versus repressed distinction. For the profile HMM, any gene with no entry in the Pfam database is given a disruption score of zero, since there is no domain template to compare against. If a mutation falls inside two overlapping domains, the higher conservation score is used as it represents the more severe potential disruption. Finally, macro F1 is used as the performance metric rather than raw accuracy. Because some subtypes may appear far more often than others, a model could achieve decent accuracy by simply predicting the most common subtype every time. Macro F1 solves this by computing a separate F1 score for each of the five cancer subtypes and averaging them equally, meaning a mistake on a rare subtype is treated as just as important as a mistake on a common subtype. This ensures the model is evaluated on how well it predicts all five subtypes, not just the most frequent ones.

## Pseudocode

### Data Loading
```{python}
load_and_align_data():

    load RNA, methylation, mutation, and clinical files

    standardize patient ID format across all four files

    keep only patients present in all four datasets

    in prototype mode: subsample to 50 patients and a small gene panel

    extract subtype labels from clinical data

    RETURN rna, methylation, mutations, labels
```
### HMM
```{python}
train_methylation_hmm(methylation):

    filter to the 10,000 most variable CpG sites
    sort sites by position along the genome

    initialize a 2-state HMM:
        state 0 = active (emission mean ~ 0.1)
        state 1 = repressed (emission mean ~ 0.8)
        transitions favor staying in the same state

    train on a subsample of 200 patients using Baum-Welch:
        E-step: estimate how likely each state is at each site
        M-step: update HMM parameters to better fit those estimates
        repeat until parameters stop changing

    if one state is almost never used, re-initialize and retrain

    RETURN trained hmm, filtered site indices


decode_regulatory_states(hmm, methylation):

    for each patient:
        run Viterbi to find the most likely sequence of hidden states
        summarize into three features:
            fraction of sites in active state
            fraction of sites in repressed state
            number of times the state switches

    RETURN table of features, one row per patient
```

### Profile HMM
```{python}
score_domain_mutations(mutations, protein_sequences, profiles):

    for each patient:

        if no mutations: assign [0, 0, 0] and move on

        for each mutation:
            align the protein sequence to its Pfam domain profile
            check if the mutation falls within the domain
            if yes: record how conserved that position is
            if no: score is zero

        summarize across all mutations:
            highest conservation score (worst disruption)
            number of mutations that hit a domain
            average conservation score

    RETURN table of features, one row per patient
```

### ML (not related to algorithms but is a personal endeavor)
```{python}
build_feature_matrix(rna, reg_features, domain_features):

    process RNA:
        log transform, standardize, keep top 500 most variable genes

    standardize regulatory and domain features to the same scale

    build four feature configurations:
        rna_only = RNA features alone
        rna_hmm = RNA + regulatory features
        rna_domain = RNA + domain features
        full = RNA + regulatory + domain features

    RETURN all four configurations

train_and_evaluate(configs, labels):

    for each feature configuration:
        for each model (logistic regression, random forest, neural network):
            run stratified 5-fold cross-validation
            record accuracy and F1 score for each fold
            average across folds

    for each model:
        compute improvement in F1 over the RNA-only baseline
        for each configuration that adds features

    RETURN results table
```

## Complexity and Bottlenecks
The clear bottleneck in this pipeline is Baum-Welch training on the methylation data. Every iteration visits every patient, every CpG site, and every pair of states. On the full unfiltered matrix this could take several hours on my laptop. Filtering down to the 10,000 most variable sites and training on a subsample of 200 patients brings the time down significantly without meaningfully changing what the model learns. Memory is a related concern since the full methylation matrix would be several GB of data, so processing patients one at a time during decoding keeps it manageable. The profile HMM scoring and classification steps are not significant bottlenecks. Breast cancer tumors have relatively few mutations, so the domain scoring runs quickly as long as the Pfam profiles are loaded once at the start rather than reloaded for every mutation. Classification on 800 patients with 500 features shouldn't take too long for logistic regression and random forest, and same for the neural network. The most important mitigation across the whole pipeline is the prototype strategy. Running everything on 50 patients and 1,000 CpG sites first means each run completes in minutes, bugs can be caught faster, and the full cohort is only attempted once the pipeline is confirmed to work correctly.

## Validation and Testing Plan
The first level of testing uses small hand-crafted examples where the correct answer is known in advance. For the HMM, a fake methylation sequence with an obvious pattern is constructed with a stretch of clearly low values, followed by a stretch of clearly high values, followed by low values again. Viterbi should decode this as active, repressed, active without ambiguity. If it doesn't, there is a bug in the decoding logic. For the profile HMM, the TP53 R175H mutation serves as a reliable test case. This is a well-studied mutation in cancer biology, and is known to be highly conserved. The model should return a high disruption score for this mutation. A score of zero would indicate the domain mapping is broken.

Once the pipeline runs correctly on toy examples it is tested on the full TCGA-BRCA cohort. Three things are checked to confirm the results are biologically sensible. First, the RNA-only model should achieve a high level of accuracy, since transcriptomics data has been stated to be a reliable method of identification. Second, patients whose HMM assigns more sites to the active state should tend to have higher overall gene expression. If this correlation is near zero, the HMM has not recovered a meaningful regulatory signal. 

A small set of automated checks will be run after each module to catch silent errors early. These include confirming that HMM transition matrix rows sum to one after each update, that patients with no mutations always receive a mutation vector of zero, that feature matrix shapes are consistent across all four configurations, and that no patient IDs are lost between modules. If any of these fail, the pipeline stops and reports where the problem occurred rather than continuing with corrupted data.

## Updated Pitfall and Risk Log
Most risks from the original proposal remain relevant and are addressed by the pipeline design. Patient overlap is handled by taking the intersection of all four datasets before any analysis begins. High dimensionality is managed by filtering to the top 500 genes and 10,000 CpG sites. The risk that RNA already captures everything is real, but is treated as a valid scientific outcome rather than a failure. Additionally, running on a MacBook M4 with 16GB RAM is a practical constraint worth acknowledging. The mitigations built into the pipeline, such as the aggressive site filtering, training on a patient subsample, and prototyping first, are designed specifically with this hardware in mind without compromising the scientific question.

## AI Apendix
Claude Sonnet 4.6 was used to understand the HMM framework and the pseudocode drafting process
