import os
import re
import pickle
import numpy as np
import pandas as pd
import pyhmmer
from Bio import Entrez, SeqIO
from io import StringIO

os.chdir(os.path.dirname(os.path.abspath(__file__)))

Entrez.email = "sshahib7@gmail.com"

#parse the maf file
def parse_maf(filepath):
    maf = pd.read_csv(filepath, sep='\t')
    print(f"Effect types present:\n{maf['effect'].value_counts()}")

    #keep only missense mutations
    maf = maf[maf['effect'] == 'Missense_Mutation'].copy()
    print(f"\nAfter missense filter: {len(maf)} mutations")

    #standardize barcodes to 12 characters
    maf['patient'] = maf['sample'].str[:12]

    #keep only patients in the cohort
    cohort_patients = pd.read_csv('Data/RNAseq.tsv', sep='\t', index_col=0).columns.str[:12].tolist()
    maf = maf[maf['patient'].isin(cohort_patients)]

    #keep only genes mutated in at least 2 patients
    gene_counts = maf.groupby('gene')['patient'].nunique()
    recurrent_genes = gene_counts[gene_counts >= 2].index
    maf = maf[maf['gene'].isin(recurrent_genes)]

    #parse protein position and polyphen score
    maf['protein_position'] = maf['Amino_Acid_Change'].apply(lambda x: int(m.group()) if (m := re.search(r'\d+', str(x))) else None)
    maf['polyphen_score'] = maf['PolyPhen'].apply(lambda x: float(m.group(1)) if (m := re.search(r'\(([\d.]+)\)', str(x))) else None)
    maf = maf.dropna(subset=['protein_position'])
    maf['protein_position'] = maf['protein_position'].astype(int)
    return maf[['patient', 'gene', 'protein_position', 'polyphen_score']]

#fetch the sequences
def fetch_sequence(gene):
    try:
        handle = Entrez.esearch(
            db="protein",
            term=f"{gene}[Gene Name] AND Homo sapiens[Organism] AND refseq[Filter]",
            retmax=1
        )
        ids = Entrez.read(handle)['IdList']
        handle.close()

        if not ids:
            return None

        handle = Entrez.efetch(
            db="protein", id=ids[0],
            rettype="fasta", retmode="text"
        )
        seq = str(SeqIO.read(StringIO(handle.read()), "fasta").seq)
        handle.close()
        return seq

    except Exception as e:
        print(f"  Warning: could not fetch sequence for {gene}: {e}")
        return None


#functions to handle cache
def build_caches(maf_df, pfam_path):
    hmm_index = pfam_path + ".h3i"
    if not os.path.exists(hmm_index):
        with pyhmmer.plan7.HMMFile(pfam_path) as hmm_file:
            pyhmmer.hmmer.hmmpress(hmm_file, pfam_path)

    with pyhmmer.plan7.HMMFile(pfam_path) as hmm_file:
        pfam_profiles = [h for h in hmm_file if h.alphabet.is_amino()]
    print(f"Loaded {len(pfam_profiles)} Pfam profiles")

    unique_genes = maf_df['gene'].unique()
    print(f"Fetching sequences for {len(unique_genes)} unique genes...")

    protein_cache = {}
    for i, gene in enumerate(unique_genes):
        protein_cache[gene] = fetch_sequence(gene)
        if (i + 1) % 20 == 0:
            print(f"  Fetched {i + 1}/{len(unique_genes)}")

    os.makedirs('Data', exist_ok=True)
    pickle.dump(protein_cache, open('Data/protein_cache.pkl', 'wb'))
    fetched = sum(v is not None for v in protein_cache.values())

    #build all the sequences in one batch for a single hmmscan call
    print("\nBuilding digital sequence batch...")
    alphabet = pyhmmer.easel.Alphabet.amino()
    digital_seqs = []
    valid_genes  = []
    for gene, seq in protein_cache.items():
        if seq is not None:
            digital_seqs.append(
                pyhmmer.easel.TextSequence(
                    name=gene.encode(),
                    sequence=seq
                ).digitize(alphabet)
            )
            valid_genes.append(gene)

    domain_cache = {gene: [] for gene in protein_cache}
    for i, hits_per_query in enumerate(
        pyhmmer.hmmscan(digital_seqs, pfam_profiles, E=1e-5, cpus=8)
    ):
        gene = valid_genes[i]
        for hit in hits_per_query:
            for domain in hit.domains:
                domain_cache[gene].append({
                    'domain': hit.name,
                    'start': domain.env_from,
                    'end': domain.env_to
                })

    pickle.dump(domain_cache, open('Data/domain_cache.pkl', 'wb'))
    hits_found = sum(len(v) > 0 for v in domain_cache.values())
    print(f"Domain cache saved — {hits_found}/{len(domain_cache)} genes had domain hits")

    return protein_cache, domain_cache

def load_caches(maf_df, pfam_path):
    if os.path.exists('Data/protein_cache.pkl') and os.path.exists('Data/domain_cache.pkl'):
        protein_cache = pickle.load(open('Data/protein_cache.pkl', 'rb'))
        domain_cache  = pickle.load(open('Data/domain_cache.pkl', 'rb'))
        print(f"  Protein cache: {len(protein_cache)} genes")
        print(f"  Domain cache:  {len(domain_cache)} genes")
        return protein_cache, domain_cache

    return build_caches(maf_df, pfam_path)



#mutation Scoring
def score_mutation(position, domain_hits, polyphen_score):
    if not domain_hits or polyphen_score is None:
        return 0.0

    scores = [polyphen_score for h in domain_hits
        if h['start'] <= position <= h['end']]
    return max(scores) if scores else 0.0


#aggregating features per patient
def compute_features(maf_df, domain_cache, all_patients):
    print(f"Computing features for {len(all_patients)} patients...")
    rows = []

    for patient in all_patients:
        mutations = maf_df[maf_df['patient'] == patient]

        if mutations.empty:
            rows.append({
                'patient': patient,
                'max_disruption': 0.0,
                'n_domain_hits': 0,
                'mean_disruption': 0.0
            })
            continue

        scores = [score_mutation(
                r.protein_position,
                domain_cache.get(r.gene, []),
                r.polyphen_score)
            for r in mutations.itertuples()]

        rows.append({
            'patient': patient,
            'max_disruption': max(scores),
            'n_domain_hits': sum(s > 0 for s in scores),
            'mean_disruption': np.mean(scores)
        })

    df = pd.DataFrame(rows).set_index('patient')
    os.makedirs('Data', exist_ok=True)
    df.to_csv('Data/domain_features.tsv', sep='\t')
    print(f"Feature matrix shape: {df.shape}")
    print(f"Saved to Data/domain_features.tsv")
    return df


if __name__ == "__main__":
    PFAM_PATH = "Data/Pfam-A.hmm"

    #parse MAF
    maf_df = parse_maf('Data/mutation.tsv')
    print(f"\nSample of parsed MAF:\n{maf_df.head()}\n")

    #load protein and domain caches
    protein_cache, domain_cache = load_caches(maf_df, PFAM_PATH)

    #get full cohort patient list from RNA file
    all_patients = pd.read_csv(
        'Data/RNAseq.tsv', sep='\t', index_col=0
    ).columns.str[:12].tolist()
    print(f"Total cohort patients: {len(all_patients)}")

    #compute per-patient domain disruption features
    features_df = compute_features(maf_df, domain_cache, all_patients)
    print(f"\nFeature preview:\n{features_df.head()}")
