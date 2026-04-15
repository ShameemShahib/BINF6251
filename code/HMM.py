import numpy as np
import pandas as pd
from scipy.stats import norm
import pickle

class HMM:
    def __init__(self):
        self.start_probs = np.array([0.5, 0.5])
        self.trans_probs = np.array([[0.9, 0.1], [0.1, 0.9]])
        self.emit_means = np.array([0.1, 0.8])
        self.emit_vars = np.array([0.01, 0.01])

    #compute emission probabilities for each state and observation
    def _emit(self, obs):
        return np.column_stack([
            norm.pdf(obs, self.emit_means[s], np.sqrt(max(self.emit_vars[s], 1e-6)))
            for s in range(2)
        ])
    #forward algorithm
    def forward(self, obs):
        E = self._emit(obs)
        T = len(obs)
        fwd = np.zeros((T, 2))
        scale = np.zeros(T)

        fwd[0] = self.start_probs * E[0]
        scale[0] = max(fwd[0].sum(), 1e-300)
        fwd[0] /= scale[0]

        for t in range(1, T):
            fwd[t] = E[t] * (fwd[t-1] @ self.trans_probs)
            scale[t] = max(fwd[t].sum(), 1e-300)
            fwd[t] /= scale[t]

        return fwd, scale, np.sum(np.log(scale))

    #backward algorithm
    def backward(self, obs, scale):
        E = self._emit(obs)
        T = len(obs)
        bwd = np.zeros((T, 2))
        bwd[-1] = 1.0 / scale[-1]

        for t in range(T-2, -1, -1):
            bwd[t] = (self.trans_probs * E[t+1] * bwd[t+1]).sum(axis=1) / max(scale[t], 1e-300)

        return bwd
    #e_step part of Baum-Welch
    def e_step(self, obs):
        fwd, scale, ll = self.forward(obs)
        bwd = self.backward(obs, scale)
        E = self._emit(obs)

        gamma = fwd * bwd
        gamma /= np.maximum(gamma.sum(axis=1, keepdims=True), 1e-300)

        xi = fwd[:-1, :, None] * self.trans_probs[None] * E[1:, None, :] * bwd[1:, None, :]
        xi /= np.maximum(xi.sum(axis=(1,2), keepdims=True), 1e-300)

        return gamma, xi, ll

    #m_step part of Baum-Welch
    def m_step(self, all_obs, all_gamma, all_xi):
        self.start_probs = np.mean([g[0] for g in all_gamma], axis=0)

        xi_sum  = sum(x.sum(axis=0) for x in all_xi)
        gamma_sum = sum(g[:-1].sum(axis=0) for g in all_gamma)
        self.trans_probs = xi_sum / np.maximum(gamma_sum[:, None], 1e-300)

        for s in range(2):
            w = np.concatenate([g[:, s] for g in all_gamma])
            o = np.concatenate(all_obs)
            den = w.sum()
            if den > 0:
                self.emit_means[s] = (w * o).sum() / den
                self.emit_vars[s]  = max((w * (o - self.emit_means[s])**2).sum() / den, 1e-4)

    #Trains the HMM using Baum-Welch
    def train(self, data, n_iter=50, tol=1e-4):
        prev_ll = -np.inf
        for it in range(n_iter):
            results = [self.e_step(obs) for obs in data]
            all_gamma, all_xi, lls = zip(*[(g, x, ll) for g, x, ll in results])
            total_ll = sum(lls)
            self.m_step(list(data), list(all_gamma), list(all_xi))
            print(f"Iter {it+1} | LL: {total_ll:.4f} | means: {self.emit_means}")
            if abs(total_ll - prev_ll) < tol:
                print("Converged")
                break
            prev_ll = total_ll

    #viterbi algorithm
    def viterbi(self, obs):
        E = np.log(self._emit(obs) + 1e-300)
        T = len(obs)
        dp = np.zeros((T, 2))
        back = np.zeros((T, 2), dtype=int)

        dp[0] = np.log(self.start_probs) + E[0]

        for t in range(1, T):
            scores = dp[t-1, :, None] + np.log(self.trans_probs + 1e-300)
            back[t] = scores.argmax(axis=0)
            dp[t]   = E[t] + scores.max(axis=0)

        states = np.zeros(T, dtype=int)
        states[-1] = dp[-1].argmax()
        for t in range(T-2, -1, -1):
            states[t] = back[t+1, states[t+1]]

        return states

    def decode_all(self, data):
        results = []
        for obs in data:
            s = self.viterbi(obs)
            results.append([np.mean(s==0), np.mean(s==1), int((s[1:]!=s[:-1]).sum())])
        return np.array(results)

    def save(self, path):
        with open(path, "wb") as f: pickle.dump(self, f)


def load(path):
        with open(path, "rb") as f: return pickle.load(f)


if __name__ == "__main__":

    #load methylation matrix
    methyl_df = pd.read_csv('Data/methyl.tsv', sep='\t', index_col=0)
    patient_ids = methyl_df.columns.tolist()

    #select top 10,000 variable sites
    site_vars = methyl_df.var(axis=1)
    top_sites = site_vars.nlargest(10000).index
    methyl_df = methyl_df.loc[top_sites]

    #align methylation data using annotation file
    annot = pd.read_csv('Data/methyl_annotation.tsv', sep='\t', index_col=0)

    common_sites = methyl_df.index.intersection(annot.index)
    methyl_df = methyl_df.loc[common_sites]
    annot = annot.loc[common_sites]

    #order sites by genomic coordinates
    chrom_order = {f'chr{i}': i for i in range(1, 23)}
    chrom_order['chrX'] = 23
    chrom_order['chrY'] = 24
    annot['chrom_num'] = annot['chrom'].map(chrom_order)

    before = len(annot)
    annot = annot.dropna(subset=['chrom_num'])
    dropped = before - len(annot)

    annot_sorted = annot.sort_values(['chrom_num', 'chromStart'])
    methyl_df = methyl_df.loc[annot_sorted.index]

    #convert to patient x site matrix and impute missing values with mean of site
    methyl_matrix = methyl_df.T.values.astype(np.float64)

    site_means = np.nanmean(methyl_matrix, axis=0)
    nan_rows, nan_cols = np.where(np.isnan(methyl_matrix))
    methyl_matrix[nan_rows, nan_cols] = site_means[nan_cols]

    #Train on randomized subset
    np.random.seed(42)
    n_patients = methyl_matrix.shape[0]
    train_idx = np.random.choice(n_patients, size=min(200, n_patients), replace=False)
    training_data = [methyl_matrix[i] for i in train_idx]

    #check the forward before going ahead will full training
    hmm_test = HMM()
    fwd, scale, ll = hmm_test.forward(training_data[0])

    #train and save the model
    hmm = HMM()
    hmm.train(training_data, n_iter=100, tol=1e-4)
    hmm.save('Data/hmm_model.pkl')

    #validate using toy sequence
    print("\nValidation test...")
    toy = np.array([0.1, 0.1, 0.1, 0.8, 0.8, 0.8])
    toy_result = hmm.viterbi(toy)
    print(f"Toy decoded: {toy_result}")
    print(f"Expected:    [0 0 0 1 1 1]")

    if not np.array_equal(toy_result, [0, 0, 0, 1, 1, 1]):
        print("validation failed")
    else:
        print("Validation passed.")

        all_data = [methyl_matrix[i] for i in range(n_patients)]
        features = hmm.decode_all(all_data)

        reg_df = pd.DataFrame(
            features,
            columns=['frac_active', 'frac_repressed', 'n_transitions'],
            index=patient_ids)
        reg_df.index.name = 'patient'
        reg_df.to_csv('Data/regulatory_features.tsv', sep='\t')

        print(f"\nSanity checks:")
        print(f"  All fractions sum to 1.0: "
              f"{np.allclose(reg_df['frac_active'] + reg_df['frac_repressed'], 1.0)}")
        print(f"  Mean frac_active:         {reg_df['frac_active'].mean():.3f}")
        print(f"  Mean frac_repressed:      {reg_df['frac_repressed'].mean():.3f}")
        print(f"  Mean n_transitions:       {reg_df['n_transitions'].mean():.1f}")
        print(f"\nSaved to Data/regulatory_features.tsv")
