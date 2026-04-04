"""
Lore KEM: Physical Decryption Failure Rate (DFR) Simulator

This script simulates the empirical Decryption Failure Rate of the Lore PKE/KEM 
by closely modeling the actual physical operations during the encryption and 
decryption processes. 

Unlike conservative theoretical bounds that assume independent errors and use 
the triangle inequality, this simulator precisely models:
1. Error Cancellation (Soft Decision): Positive and negative polynomial noise 
   canceling out in Repetition encoding aggregation.
2. Anti-cyclic Convolution Sign Flipping: Modeled via 50/50 mixed distributions.
3. Burst Error Correlation: Integrated full-probability models conditioning on 
   the base noise energy (non-zero terms) to capture highly correlated polynomial 
   coefficient failures under BCH encoding.

Dependencies:
    - proba_util.py (Requires `iter_law_convolution`, `law_convolution`)
"""

import math
import sys

from proba_util import (
    iter_law_convolution,
    law_convolution,
    tail_probability
)

# ==============================================================================
# Statistical Utility Functions
# ==============================================================================
def binomial_tail_prob(n_trials, p_success, threshold):
    """Computes the upper tail probability P(X > threshold) for a Binomial distribution."""
    prob = 0.0
    for k in range(threshold + 1, n_trials + 1):
        try:
            comb = math.comb(n_trials, k)
        except AttributeError:
            comb = math.factorial(n_trials) // (math.factorial(k) * math.factorial(n_trials - k))
        
        if comb == float('inf'): 
            break
            
        term = comb * (p_success ** k) * ((1 - p_success) ** (n_trials - k))
        prob += term
        
        if term < 1e-100 and prob > 0: 
            break
            
    return prob

def binomial_pmf(n_trials, k, p):
    """Computes the Binomial PMF P(X=k) using log-gamma to prevent overflow."""
    if p == 0: 
        return 1.0 if k == 0 else 0.0
    if p == 1: 
        return 1.0 if k == n_trials else 0.0
    if k < 0 or k > n_trials: 
        return 0.0
        
    log_prob = (math.lgamma(n_trials + 1) - math.lgamma(k + 1) - math.lgamma(n_trials - k + 1) 
                + k * math.log(p) + (n_trials - k) * math.log(1 - p))
    return math.exp(log_prob)

def condition_law(base_law, rho):
    """
    Reshapes the error distribution based on a given non-zero probability (rho).
    Crucial for simulating base noise energy fluctuations to compute full-probability 
    integrals for burst error correlations in BCH codes.
    """
    p_nonzero_orig = sum(prob for val, prob in base_law.items() if val != 0)
    if p_nonzero_orig == 0: 
        return base_law

    cond_law = {}
    for val, prob in base_law.items():
        if val == 0:
            cond_law[val] = 1.0 - rho
        else:
            cond_law[val] = prob * (rho / p_nonzero_orig)
            
    return cond_law

def get_best_bch_config(available_space, plaintext_bits):
    """Determines the optimal BCH configuration based on available polynomial space."""
    m_down = 0
    while (2 ** (m_down + 1) - 1) <= available_space: 
        m_down += 1
        
    t_down = -1
    n_trials_down = 0
    if m_down > 0:
        n_std_down = 2 ** m_down - 1
        if plaintext_bits <= n_std_down:
            parity_down = n_std_down - plaintext_bits
            t_down = parity_down // m_down
            n_trials_down = n_std_down
            
    m_up = 1
    while (2 ** m_up - 1) < available_space: 
        m_up += 1
        
    parity_up = available_space - plaintext_bits
    t_up = parity_up // m_up if parity_up >= 0 else -1
    n_trials_up = available_space

    if t_down > t_up:
        return t_down, m_down, n_trials_down, "Truncated (Downwards)"
    return t_up, m_up, n_trials_up, "Shortened (Upwards)"

# ==============================================================================
# Base Noise Distribution Builders
# ==============================================================================
def build_cu_law_pure_2i_plus_1(t):
    """Builds the CRT compression error distribution for public key and Cu."""
    error_counts = {}
    S_R = [2 * i + 1 for i in range(t) if 2 * i + 1 < t]
    
    for curr_t in range(t):
        candidates = []
        for r in S_R:
            diff = (r - curr_t) % t
            if diff > t // 2: 
                diff -= t
                
            if t % 2 == 0 and abs(diff) * 2 == t:
                candidates.extend([{'diff': abs(diff), 'abs': abs(diff)}, 
                                   {'diff': -abs(diff), 'abs': abs(diff)}])
            else:
                candidates.append({'diff': diff, 'abs': abs(diff)})
                
        min_abs = min(c['abs'] for c in candidates)
        best_valid = sorted([c for c in candidates if c['abs'] == min_abs], key=lambda x: -x['diff'])
        err = best_valid[0]['diff']
        error_counts[err] = error_counts.get(err, 0) + 1.0
        
    total = sum(error_counts.values())
    return {k: v / total for k, v in error_counts.items()}

def build_cv_law_pure_compress(q, bits=1):
    """Builds the CRT compression error distribution for Cv."""
    error_counts = {}
    num_intervals = 2 ** bits
    for xq in range(q):
        interval_idx = round(xq * num_intervals / q) % num_intervals
        center_q = round(interval_idx * q / num_intervals)
        diff = center_q - xq
        
        if diff > q // 2: diff -= q
        elif diff < -q // 2: diff += q
            
        error_counts[diff] = error_counts.get(diff, 0) + 1.0
        
    total = sum(error_counts.values())
    return {k: v / total for k, v in error_counts.items()}

def build_e_prime_prime_law(t):
    """Builds the uniform noise distribution for e''."""
    return {v: 1.0 / t for v in range(-t // 2 + 1, t // 2 + 1)}

# ==============================================================================
# Core Simulation Pipeline
# ==============================================================================
def compute_pipeline(law_cu_current, config):
    """
    Computes the final error distribution for a single coefficient.
    Features a 50/50 mixture to simulate sign flipping caused by the X^n = -1 
    anti-cyclic polynomial multiplication over the ring.
    """
    # 1. Prepare mirrored and scaled noise distributions
    law_cu_neg = {-v: p for v, p in law_cu_current.items()}
    law_cu_2x = {2 * v: p for v, p in law_cu_current.items()}
    law_cu_neg_2x = {-2 * v: p for v, p in law_cu_current.items()}

    # 2. 50/50 Mixture modeling for mean-drift stabilization
    law_cu_mixed_1 = {}
    for v, p in law_cu_current.items(): law_cu_mixed_1[v] = law_cu_mixed_1.get(v, 0) + p * 0.5
    for v, p in law_cu_neg.items(): law_cu_mixed_1[v] = law_cu_mixed_1.get(v, 0) + p * 0.5

    law_cu_mixed_2 = {}
    for v, p in law_cu_2x.items(): law_cu_mixed_2[v] = law_cu_mixed_2.get(v, 0) + p * 0.5
    for v, p in law_cu_neg_2x.items(): law_cu_mixed_2[v] = law_cu_mixed_2.get(v, 0) + p * 0.5

    # 3. Global allocation: total items spanning all k polynomials
    k = config['k']
    total_ones = (config['num_ones_per_poly'] + config['num_minus_ones_per_poly']) * k
    total_twos = (config['num_twos_per_poly'] + config['num_minus_twos_per_poly']) * k

    # 4. Strict iterative convolution to prevent false variance
    law_part_1 = iter_law_convolution(law_cu_mixed_1, total_ones)
    law_part_2 = iter_law_convolution(law_cu_mixed_2, total_twos)
    law_eT_S_prime = law_convolution(law_part_1, law_part_2)

    # 5. Convolution of the opposing secret noise (-e2^T s)
    law_sT_E_prime_neg = {-val: prob for val, prob in law_eT_S_prime.items()}
    law_intermediate = law_convolution(law_eT_S_prime, law_sT_E_prime_neg)

    # 6. Include explicit noise e'' and compression error Cv
    law_e_pp = build_e_prime_prime_law(config['t'])
    law_cv_final = build_cv_law_pure_compress(config['q'], bits=config['cv_bits'])
    
    law_intermediate = law_convolution(law_intermediate, law_e_pp)
    law_final_current = law_convolution(law_intermediate, law_cv_final)

    threshold = (config['t'] * config['q']) / 4.0
    raw_fail = sum(prob for val, prob in law_final_current.items() if abs(val) >= threshold)
    
    return law_final_current, raw_fail

# ==============================================================================
# Main Execution
# ==============================================================================
if __name__ == "__main__":
    # Scheme Configuration (Currently aligned with User's target parameters)
    CONFIG = {
        'n': 768,
        'k': 3,
        'q': 257,
        't': 4,
        'cv_bits': 1,
        'plaintext_len': 512,
        'num_twos_per_poly': 20,
        'num_minus_twos_per_poly': 20,
        'num_ones_per_poly': 140,
        'num_minus_ones_per_poly': 140
    }
    
    AVAILABLE_SPACE = CONFIG['n']
    THRESHOLD = (CONFIG['t'] * CONFIG['q']) / 4.0

    print("--- Lore Scheme Parameters (Physical Simulation Model) ---")
    
    if AVAILABLE_SPACE % CONFIG['plaintext_len'] == 0 and (AVAILABLE_SPACE / CONFIG['plaintext_len']) >= 2:
        K_VEC = AVAILABLE_SPACE // CONFIG['plaintext_len']
        USE_ECC = False
        print(f"[Strategy] Space ({AVAILABLE_SPACE}) is multiple of Msg ({CONFIG['plaintext_len']}).")
        print(f"           -> Using K_VEC mode. k_vec = {K_VEC}")
    else:
        K_VEC = 1
        USE_ECC = True
        bch_t, bch_m, n_trials, strategy_name = get_best_bch_config(AVAILABLE_SPACE, CONFIG['plaintext_len'])
        print(f"[Strategy] Space ({AVAILABLE_SPACE}) vs Msg ({CONFIG['plaintext_len']}).")
        print(f"           -> Using Best BCH Strategy: {strategy_name}")
        print(f"           -> Effective Coding Length (n_trials): {n_trials}")
        print(f"           -> Galois Field: m={bch_m}")
        print(f"           -> Correctable Errors (t) = {bch_t}")

    print("-" * 60 + "\n")

    # 1. Base DFR Calculation
    law_cu_base = build_cu_law_pure_2i_plus_1(CONFIG['t'])
    law_final_base, raw_fail_prob = compute_pipeline(law_cu_base, CONFIG)

    print("--- Final DFR Analysis ---")
    print(f"Base Threshold: {THRESHOLD}")
    print(f"Raw Single Coefficient Fail Prob (p): {raw_fail_prob}")
    if raw_fail_prob > 0:
        print(f"   ≈ 2^({math.log2(raw_fail_prob):.2f})")

    # 2. Error Correction & Aggregation Evaluation
    if not USE_ECC:
        print(f"\n[Mode: K_VEC = {K_VEC}]")
        if K_VEC == 1:
            final_dfr = CONFIG['plaintext_len'] * raw_fail_prob
        else:
            # Applies physical error cancellation before absolute evaluation
            law_vec = iter_law_convolution(law_final_base, K_VEC)
            theta_vec = K_VEC * THRESHOLD
            vec_fail_prob = sum(prob for val, prob in law_vec.items() if abs(val) >= theta_vec)

            print(f"Vector (Block) Fail Prob: {vec_fail_prob}")
            final_dfr = CONFIG['plaintext_len'] * vec_fail_prob

        final_dfr = min(1.0, final_dfr)
        print(f"Final Scheme DFR: {final_dfr}")
        if final_dfr > 0:
            print(f"   ≈ 2^({math.log2(final_dfr):.2f})")

    else:
        print(f"\n[Mode: Optimized BCH - {strategy_name}]")
        print(f"Code Params: Trials(N)={n_trials}, Msg(K)={CONFIG['plaintext_len']}, T={bch_t}")

        if raw_fail_prob == 0:
            ecc_fail_prob = 0.0
        else:
            # Full-probability integration for burst error correlation
            N_noise = 2 * CONFIG['n'] * CONFIG['k']
            p_nonzero_base = sum(prob for val, prob in law_cu_base.items() if val != 0)

            if p_nonzero_base == 0.0 or p_nonzero_base == 1.0:
                ecc_fail_prob = binomial_tail_prob(n_trials, raw_fail_prob, bch_t)
            else:
                mu_W = N_noise * p_nonzero_base
                sigma_W = math.sqrt(N_noise * p_nonzero_base * (1 - p_nonzero_base))

                final_bch_dfr = 0.0
                min_W = max(0, int(mu_W - 6 * sigma_W))
                max_W = min(N_noise, int(mu_W + 6 * sigma_W))
                step = max(1, int(sigma_W / 10))

                for w in range(min_W, max_W + 1, step):
                    p_w = binomial_pmf(N_noise, w, p_nonzero_base)
                    if step > 1: 
                        p_w *= step
                    if p_w < 1e-100: 
                        continue

                    # Condition the noise on the current burst state
                    rho = w / N_noise
                    law_cu_cond = condition_law(law_cu_base, rho)
                    _, p_fail_cond = compute_pipeline(law_cu_cond, CONFIG)

                    # Accumulate weighted probabilities
                    bch_fail_cond = binomial_tail_prob(n_trials, p_fail_cond, bch_t)
                    final_bch_dfr += p_w * bch_fail_cond

                ecc_fail_prob = final_bch_dfr

        print(f"ECC Failure Prob (Integrated over Burst Errors): {ecc_fail_prob}")
        if ecc_fail_prob > 0:
            print(f"   ≈ 2^({math.log2(ecc_fail_prob):.2f})")

    print("\nSimulation Complete.")