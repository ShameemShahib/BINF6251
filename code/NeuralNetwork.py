import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import f1_score, accuracy_score


#loading features from previous scripts
def load_features():
    rna = pd.read_csv('Data/RNAseq.tsv',             sep='\t', index_col=0).T
    reg = pd.read_csv('Data/regulatory_features.tsv', sep='\t', index_col=0)
    dom = pd.read_csv('Data/domain_features.tsv',     sep='\t', index_col=0)
    clin = pd.read_csv('Data/clinicalMatrix.tsv',      sep='\t', index_col=0)
    labels = clin['PAM50Call_RNAseq'].dropna()

    #standardize all indices to 12 characters
    rna.index = rna.index.str[:12]
    reg.index = reg.index.str[:12]
    dom.index = dom.index.str[:12]
    labels.index = labels.index.str[:12]

    #drop any duplicates
    rna = rna[~rna.index.duplicated(keep='first')]
    reg = reg[~reg.index.duplicated(keep='first')]
    dom = dom[~dom.index.duplicated(keep='first')]
    labels = labels[~labels.index.duplicated(keep='first')]

    #align to shared patients
    patients = sorted(
        set(rna.index) & set(reg.index) & set(dom.index) & set(labels.index))
    print(f"Shared patients after alignment: {len(patients)}")

    rna, reg, dom, labels = (
        rna.loc[patients],
        reg.loc[patients],
        dom.loc[patients],
        labels.loc[patients])

    print(f"Subtype counts:\n{labels.value_counts()}")

    #preprocess RNA
    rna = np.log1p(rna)
    rna = rna[rna.var(axis=0).nlargest(500).index]
    rna = (rna - rna.mean()) / rna.std()
    rna = rna.fillna(0)

    #standardize regulatory and domain features
    reg = (reg - reg.mean()) / reg.std()
    dom = (dom - dom.mean()) / dom.std()

    #fill any missing values
    reg = reg.fillna(0).replace([np.inf, -np.inf], 0)
    dom = dom.fillna(0).replace([np.inf, -np.inf], 0)

    print(f"Domain feature distributions:")
    print(dom.describe())

    configs = {
        'rna_only':   rna,
        'rna_hmm':    pd.concat([rna, reg], axis=1),
        'rna_domain': pd.concat([rna, dom], axis=1),
        'full':       pd.concat([rna, reg, dom], axis=1)
    }

    return configs, labels

configs, labels = load_features()


#Network
class SubtypeClassifier(nn.Module):
    def __init__(self, n_features):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_features, 64),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(64, 5)
        )
    def forward(self, x):
        return self.net(x)

#training Loop
def run_fold(X_train, y_train, X_val, y_val, n_features, n_epochs=100):
    model = SubtypeClassifier(n_features)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    criterion = nn.CrossEntropyLoss()

    model.train()
    for _ in range(n_epochs):
        optimizer.zero_grad()
        loss = criterion(model(torch.FloatTensor(X_train)), torch.LongTensor(y_train))
        loss.backward()
        optimizer.step()

    model.eval()
    with torch.no_grad():
        preds = model(torch.FloatTensor(X_val)).argmax(dim=1).numpy()

    return f1_score(y_val, preds, average='macro'), accuracy_score(y_val, preds)

#cross-Validation
def evaluate(configs, labels):
    label_map = {'LumA': 0, 'LumB': 1, 'Her2': 2, 'Basal': 3, 'Normal': 4}

    #map labels and drop any that didn't match
    y_series = labels.map(label_map)
    unrecognized = labels[y_series.isna()].unique()

    valid_mask = y_series.notna()
    y = y_series[valid_mask].values.astype(int)

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    results = {}

    for name, features in configs.items():
        X = features.values[valid_mask][:]
        X = X.astype(np.float32)

        fold_scores = [
            run_fold(X[tr], y[tr], X[va], y[va], X.shape[1])
            for tr, va in skf.split(X, y)]

        f1s, accs = zip(*fold_scores)
        results[name] = {
            'mean_f1': np.mean(f1s),
            'std_f1': np.std(f1s),
            'mean_acc': np.mean(accs)}
        print(f"{name}: F1={np.mean(f1s):.3f} ± {np.std(f1s):.3f} | Acc={np.mean(accs):.3f}")

    return results

#summarize results
def summarize(results):
    baseline = results['rna_only']['mean_f1']
    print("\n=== Improvement over RNA-only baseline ===")
    for name, m in results.items():
        print(f"{name}: F1={m['mean_f1']:.3f} (±{m['std_f1']:.3f}) | "
              f"improvement={m['mean_f1'] - baseline:+.3f}")

if __name__ == "__main__":
    configs, labels = load_features()
    results = evaluate(configs, labels)
    summarize(results)
