# Integrating DNA Methylation and Protein Domain Features to Improve Breast Cancer Subtype Prediction

## Research Question
#### Background
Gene expression profiling is one of the most reliable ways to identify tumor subtype. In breast cancer, transcriptomic patterns clearly separate tumors into biologically meaningful groups. However, RNA mainly tells us what the cell is doing at a given moment. It does not directly explain how genes are being regulated or whether key proteins are functionally damaged. In other words, transcriptomic data provide a detailed snapshot, but not the full story.

This project investigates whether adding biological context improves subtype prediction. Specifically, I will incorporate two additional layers: regulatory information inferred from DNA methylation and functional information about whether conserved protein domains are disrupted by mutations. Rather than feeding raw data directly into a model, I will first summarize them into biologically meaningful features, such as likely active or repressed regions and potential domain disruption.

Multi-omics approaches are increasingly common in cancer research. Many methods combine raw data types and allow models to learn patterns automatically. While effective, this can make it difficult to understand what each data type contributes.

#### Question:

**If we give the model information about gene regulation and protein function, does it predict the cancer subtype better than using RNA data alone?**

By comparing models that exclusively predict based on the transcriptome with models that include this additional data, I aim to understand whether the extra layers offer new information or whether their effects are already reflected in the transcriptome.

---
## Algorithm an Algorithm Class
This project uses Hidden Markov Models (HMMs) and profile Hidden Markov Models to convert raw data into interpretable summaries. I will also compare traditional predictive models with a small neural network.

#### HMM
I will use an HMM to infer hidden regulatory states from DNA methylation data. Methylation values are noisy measurements, while regulatory activity (such as “active” or “repressed”) is not directly observable. HMMs are well suited for this problem because they infer hidden states that explain observed data. The output will be simplified summaries, such as the proportion of regions likely active or repressed in each tumor.

#### Profile HMM
I will use profile HMMs to evaluate whether mutations affect conserved protein domains. A mutation alone does not indicate functional impact. Profile HMMs allow me to assess how well a protein sequence matches known domain patterns, helping estimate whether key functional regions are likely disrupted.

#### Neural Network (Experimental)
I will also build a small neural network that takes RNA, regulatory summaries, and domain features as input. Each data type will first be processed separately and then combined for subtype prediction. Neural networks can capture complex relationships, allowing me to test whether the interactions between these biological layers influence performance. The architecture will remain simple to keep the focus on feature evaluation rather than model optimization.

---
## Data Plan
#### Source
Data will come from The Cancer Genome Atlas (TCGA), which provides well-curated multi-omics datasets for breast cancer. TCGA includes subtype labels and multiple molecular measurements across overlapping patients, making it well-suited for this project.

#### Data Types and Acquisition
For each patient, I will collect gene expression, DNA methylation, somatic mutation calls, and subtype labels. Gene expression will serve as the baseline. Methylation data will be used in the HMM to infer regulatory states. Mutation data will be mapped to protein sequences and evaluated for domain disruption using profile HMMs. All datasets are publicly available and de-identified, so no special access is required.

#### Prototype Data Plan
I will begin with a smaller subset of patients and a reduced gene list to test the pipeline. This stage will confirm that data alignment, HMM inference, and domain mapping work correctly. Once validated, the same framework will be scaled to a larger dataset. This staged approach keeps the project manageable while preserving the overall design.

---
## Success Criteria
Success means building a reproducible pipeline and clearly evaluating whether regulatory and protein domain features add value beyond RNA alone. The goal is not just higher accuracy, but understanding whether additional biological layers provide complementary information. Success in this case is a model that can show how influential an omics dataset is, or if it has no influence on prediction at all.

I will compare RNA-only models with models that include added features and report metrics such as accuracy and F1 score. I will also examine whether certain subtypes show stronger associations with regulatory or functional features.

To ensure results are reasonable, RNA-only models should perform well, since expression is known to distinguish subtypes. I will also check for expected biological patterns, such as agreement between inferred regulatory activity and gene expression. If results contradict known biology, preprocessing or modeling steps will be revisited. Consistent trends across simple models and neural networks will increase confidence in the findings.

---
## Pitfall Scan
Multi-omics analysis presents several challenges. One concern is incomplete data overlap. Not all patients have all omics types, which may reduce sample size. I will assess overlap early and define inclusion criteria accordingly. High dimensionality is another issue, especially for RNA data, where genes outnumber patients. I will reduce dimensionality by filtering genes based on variance or biological relevance and use cross-validation to monitor overfitting. HMM states may not always align neatly with intuitive biological categories. I will compare inferred regulatory states with gene expression patterns and adjust the model if inconsistencies arise. Computational limitations are also possible, particularly for protein domain searches on a laptop (a MacBook in my case). I will start with curated gene sets and expand only if the runtime is manageable. Finally, RNA may already capture much of the subtype signal. If additional layers produce only small improvements, this will be interpreted as evidence of redundancy rather than failure. Overfitting will be controlled through cross-validation and simple model design.

---
## Planned Repository Structure
To maintain organization and reproducibility, the project will follow a modular structure separating data processing, algorithms, feature construction, and modeling.

I anticipate organizing the repository in the following way:

project/
  
  data/
  
    raw/ – original downloaded TCGA files
    processed/ – cleaned and aligned datasets

  algorithms/
  
    hmm/ – regulatory state inference from methylation
    profile/ – mutation-to-domain mapping
    features/ – construction of patient-level feature tables
    models/ – training and evaluation scripts

  results/ – saved model outputs and evaluation metrics
  
  figures/ – plots/graphs
  
  docs/ – README, proposal, and other documentation
