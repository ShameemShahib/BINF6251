import numpy as np
import sys
from HMM import HMM

#Test 1: Toy sequence decoding to ensure viterbi backtracking logic works
def test_viterbi_toy_sequence():
    hmm = HMM()
    obs = np.array([0.05, 0.08, 0.10, 0.85, 0.90, 0.88])
    states = hmm.viterbi(obs)
    assert list(states) == [0, 0, 0, 1, 1, 1], (
        f"Expected [0,0,0,1,1,1], got {list(states)}")

#Test 2: Forward pass to ensure no NaN log-likelihood
def test_forward():
    hmm = HMM()
    np.random.seed(0)
    obs = np.random.beta(0.5, 0.5, size=100)
    _, _, ll = hmm.forward(obs)
    assert np.isfinite(ll), f"Log-likelihood should be finite, got {ll}"


#Test 3: Ensure forward probabilities sum to 1 at every site
def test_forward_prob_sum():
    hmm = HMM()
    obs = np.array([0.1, 0.2, 0.8, 0.9, 0.1])
    fwd, _, _ = hmm.forward(obs)
    row_sums = fwd.sum(axis=1)
    assert np.allclose(row_sums, 1.0, atol=1e-6), (
        f"Forward prob rows should sum to 1, got {row_sums}")


#Test 4: E-step validation, check if gamma sums to 1
def test_e_step():
    hmm = HMM()
    obs = np.array([0.1, 0.5, 0.8, 0.2, 0.9])
    gamma, _, _ = hmm.e_step(obs)
    row_sums = gamma.sum(axis=1)
    assert np.allclose(row_sums, 1.0, atol=1e-6), (
        f"Gamma rows must sum to 1, got {row_sums}")


#Test 5: Ensure transition matrix rows sum to 1 after M-step
def test_transition_matrix():
    hmm = HMM()
    np.random.seed(42)
    data = [
        np.array([0.1, 0.1, 0.8, 0.8, 0.1]),
        np.array([0.8, 0.8, 0.1, 0.1, 0.8])]
    results   = [hmm.e_step(obs) for obs in data]
    all_gamma = [r[0] for r in results]
    all_xi    = [r[1] for r in results]
    hmm.m_step(data, all_gamma, all_xi)
    row_sums = hmm.trans_probs.sum(axis=1)
    assert np.allclose(row_sums, 1.0, atol=1e-6), (
        f"Transition matrix rows must sum to 1, got {row_sums}")


#Test 6: Ensure log-likelihood increases during training
def test_log_likelihood():
    hmm = HMM()
    np.random.seed(0)
    data = []
    for _ in range(10):
        active = np.random.normal(0.1, 0.05, 50).clip(0, 1)
        repressed = np.random.normal(0.8, 0.05, 50).clip(0, 1)
        obs = np.concatenate([active, repressed, active])
        data.append(obs)
    log_likelihoods = []
    for _ in range(5):
        results = [hmm.e_step(obs) for obs in data]
        all_gamma = [r[0] for r in results]
        all_xi = [r[1] for r in results]
        total_ll = sum(r[2] for r in results)
        log_likelihoods.append(total_ll)
        hmm.m_step(data, all_gamma, all_xi)

    for i in range(1, len(log_likelihoods)):
        assert log_likelihoods[i] >= log_likelihoods[i-1] - 1e-3, (
            f"Log-likelihood decreased at iteration {i+1}: "
            f"{log_likelihoods[i-1]:.4f} → {log_likelihoods[i]:.4f}")

# Test 7: Ensure there are no NaN in the data
def test_NaN():
    matrix = np.array([
        [0.1, np.nan, 0.8],
        [np.nan, 0.5, 0.9],
        [0.2, 0.3, np.nan]
    ])

    site_means = np.nanmean(matrix, axis=0)
    nan_rows, nan_cols = np.where(np.isnan(matrix))
    matrix[nan_rows, nan_cols] = site_means[nan_cols]

    assert not np.isnan(matrix).any(), "NaNs remain after imputation"
    assert np.isclose(matrix[0, 1], np.mean([0.5, 0.3])), (
        "Imputed value should equal the site mean across non-NaN patients")



if __name__ == "__main__":

    tests = [
        test_viterbi_toy_sequence,
        test_forward,
        test_forward_prob_sum,
        test_e_step,
        test_transition_matrix,
        test_log_likelihood,
        test_NaN,]

    passed = 0
    failed = 0

    print("Running HMM tests")

    for test in tests:
        name = test.__name__
        test()
        print(f"  PASSED  {name}")
        passed += 1

    print(f"Results: {passed} passed, {failed} failed out of {len(tests)} tests")
