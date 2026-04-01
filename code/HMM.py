import numpy as np
import pandas as pd
from scipy.stats import norm

def baum_welch(sequences):
    T_mat = np.array([[0.9, 0.1], [0.1, 0.9]])
    means = np.array([0.1, 0.8])

    for _ in range(20):
        T_acc = np.zeros((2, 2))
        mean_num = np.zeros(2)
        mean_den = np.zeros(2)

        for obs in sequences:

            # forward pass
            forward_prob = np.zeros((len(obs), 2))
            forward_prob[0] = 0.5 * norm.pdf(obs[0], means, 0.1)
            for t in range(1, len(obs)):
                forward_prob[t] = (forward_prob[t-1] @ T_mat) * norm.pdf(obs[t], means, 0.1)
                forward_prob[t] /= forward_prob[t].sum()

            # backward pass
            backward_prob = np.ones((len(obs), 2))
            for t in range(len(obs) - 2, -1, -1):
                backward_prob[t] = T_mat @ (norm.pdf(obs[t+1], means, 0.1) * backward_prob[t+1])
                backward_prob[t] /= backward_prob[t].sum()

            # state probabilities
            state_probs = forward_prob * backward_prob
            state_probs = state_probs / state_probs.sum(axis=1, keepdims=True)

            # accumulate
            T_acc += forward_prob[:-1].T @ (backward_prob[1:] * norm.pdf(obs[1:, None], means, 0.1))
            mean_num += (state_probs * obs[:, None]).sum(axis=0)
            mean_den += state_probs.sum(axis=0)

        # update parameters
        T_mat = T_acc / T_acc.sum(axis=1, keepdims=True)
        means = mean_num / mean_den

    return T_mat, means


def viterbi(obs, T_mat, means):
    best_scores = np.log(norm.pdf(obs[0], means, 0.1))
    log_T = np.log(T_mat)
    traceback  = np.zeros((len(obs), 2), dtype=int)

    for t in range(1, len(obs)):
        candidate_scores = best_scores + log_T.T
        traceback[t] = candidate_scores.argmax(axis=1)
        best_scores = candidate_scores.max(axis=1) + np.log(norm.pdf(obs[t], means, 0.1))

    state_sequence = np.zeros(len(obs), dtype=int)
    state_sequence[-1] = best_scores.argmax()
    for t in range(len(obs) - 2, -1, -1):
        state_sequence[t] = traceback[t+1, state_sequence[t+1]]

    return state_sequence


if __name__ == '__main__':

    # load methylation data
    methyl = pd.read_csv('Data/methyl.tsv', sep='\t', index_col=0)

    # filter to top 1000 most variable CpG sites
    top_sites = methyl.var(axis=1).nlargest(1000).index
    methyl = methyl.loc[top_sites]

    # prototype: use only 50 patients
    methyl = methyl.iloc[:, :50]

    # build training sequences from first 20 patients
    sequences = []
    for i in range(20):
        obs = methyl.iloc[:, i].dropna().values.astype(float)
        sequences.append(obs)

    # train HMM
    T_mat, means = baum_welch(sequences)

    # decode all patients and record fraction of active sites
    results = {}
    for patient in methyl.columns:
        obs = methyl[patient].dropna().values.astype(float)
        path = viterbi(obs, T_mat, means)
        results[patient] = (path == 0).mean()

    # save results
    reg_features = pd.Series(results, name='frac_active')
    reg_features.to_csv('Data/regulatory_features.tsv', sep='\t')
    print(reg_features.describe())
