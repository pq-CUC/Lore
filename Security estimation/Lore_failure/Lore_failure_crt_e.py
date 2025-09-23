
import math

# Import all necessary probability calculation functions from the proba_util.py file
from proba_util import (
    build_centered_binomial_law,
    iter_law_convolution,
    law_convolution,
    tail_probability
)


# ==============================================================================
# Helper functions (unchanged)
# ==============================================================================
def build_reconciliation_law(t):
    """
    (Original version) Build a non-uniform error distribution based on the t value.
    Used for e and E'.
    """
    if t == 3:
        return {-1: 1 / 3, 0: 1 / 3, 1: 1 / 3}
    elif t == 7:
        return {-2: 1 / 14, -1: 2 / 7, 0: 2 / 7, 1: 2 / 7, 2: 1 / 14}
    elif t == 13:
        return {-2: 1 / 26, -1: 4 / 13, 0: 4 / 13, 1: 4 / 13, 2: 1 / 26}
    else:
        raise ValueError("Unsupported t value for reconciliation. Please choose 3, 7, or 13.")


def build_uniform_law(min_val, max_val):
    """
    Build a discrete uniform distribution over the integer range [min_val, max_val].
    Used for the base distribution of e''.
    """
    num_points = max_val - min_val + 1
    prob = 1.0 / num_points
    return {i: prob for i in range(min_val, max_val + 1)}


# ==============================================================================
# Scheme parameter definitions
# ==============================================================================
print("--- Scheme Parameters ---")
n = 256  # Degree of the polynomial ring
k = 2  # Dimension of vectors/matrices (module rank)
q = 257  # Main modulus
eta = 1  # Base noise distribution parameter (only effective when RECONCILIATION_T = 0)
t_repetition = 3  # Number of scheme repetitions, only affects the decryption threshold

# Set to 3, 7, or 13 to use the new reconciliation error distribution model
# Set to 0 to use the original centered binomial distribution model
RECONCILIATION_T = 3

# The source of e, E', e'' depends on RECONCILIATION_T
i1 = 1  # Aggregation count for e (only effective when RECONCILIATION_T = 0)
i2 = 1  # Aggregation count for E' (only effective when RECONCILIATION_T = 0)
i3 = 1  # Aggregation count for e'' (only effective when RECONCILIATION_T = 0)

# Define new FHW distribution parameters
# Distribution: {-2: 1/16, -1: 1/4, 0: 3/8, 1: 1/4, 2: 1/16}
num_twos = 0
num_ones = 40
num_zeros = 176
num_minus_ones = 40
num_minus_twos = 0

print(f"n = {n}, k = {k}, q = {q}")
print(f"Repetition count t_repetition = {t_repetition}")
if RECONCILIATION_T in [3, 7, 13]:
    print(f"All noise (e, E') comes from the original (non-uniform) reconciliation error distribution for t={RECONCILIATION_T}")
    # === Start of change: Update print info for e'' ===
    print(f"Noise (e'') comes from two independent additions of the uniform distribution corresponding to t={RECONCILIATION_T}")
    # === End of change ===
else:
    print(f"All noise (e, E', e'') comes from the centered binomial distribution (η = {eta})")
    print(f"Aggregation count for e: i1 = {i1}")
    print(f"Aggregation count for E': i2 = {i2}")
    print(f"Aggregation count for e'': i3 = {i3}")
print(f"s, S' come from the new 5-value FHW distribution (±2:{num_twos}, ±1:{num_ones}, 0:{num_zeros})")
print("-" * 20 + "\n")

# ==============================================================================
# Start calculating the probability distribution of the total noise E_total = e^T*S' - s^T*E' + e''
# ==============================================================================

# 1. Generate base noise distributions
print("Step 1: Generate base noise distributions...")
if RECONCILIATION_T in [3, 7, 13]:
    # --- e and E' use the original, non-uniform distribution ---
    print(f"   -> e, E' use the original reconciliation error model for t={RECONCILIATION_T}.")
    reconciliation_law = build_reconciliation_law(RECONCILIATION_T)
    law_e = reconciliation_law
    law_E_prime = reconciliation_law

    # === Start of change: Implement the "two independent uniform additions" logic for e'' ===
    print(f"   -> e'' uses the uniform distribution model for t={RECONCILIATION_T}, summed from two independent samples.")
    # 1. Determine the bounds of the uniform distribution based on t
    if RECONCILIATION_T == 3:
        bound = 1  # Range [-1, 1]
    elif RECONCILIATION_T == 7:
        bound = 3  # Range [-3, 3]
    else:  # RECONCILIATION_T == 13
        bound = 6  # Range [-6, 6]

    # 2. Create the base uniform distribution
    base_uniform_law = build_uniform_law(-bound, bound)

    # 3. Convolve the base uniform distribution twice to get the final distribution for e''
    law_e_double_prime = iter_law_convolution(base_uniform_law, 2)
    # === End of change ===

else:
    print(f"   -> Using centered binomial distribution model (η = {eta}).")
    law_chi = build_centered_binomial_law(eta)
    law_e = iter_law_convolution(law_chi, i1)
    law_E_prime = iter_law_convolution(law_chi, i2)
    law_e_double_prime = iter_law_convolution(law_chi, i3)

print("Distributions for e, E', e'' calculated.\n")

# 2. Calculate the distributions for the two MLWE product terms (using the new 5-value FHW addition model)
print("Step 2: Calculate distributions for the MLWE product terms separately...")

# --- Calculate the distribution of e^T*S' ---
print("  --- 2a: Calculating distribution of e^T*S' ---")
law_e_neg = {-val: prob for val, prob in law_e.items()}
law_e_2x = {2 * val: prob for val, prob in law_e.items()}
law_e_neg_2x = {-2 * val: prob for val, prob in law_e.items()}
law_sum_pos_2_S = iter_law_convolution(law_e_2x, num_twos)
law_sum_pos_1_S = iter_law_convolution(law_e, num_ones)
law_sum_neg_1_S = iter_law_convolution(law_e_neg, num_minus_ones)
law_sum_neg_2_S = iter_law_convolution(law_e_neg_2x, num_minus_twos)
law_poly_product_S_1 = law_convolution(law_sum_pos_2_S, law_sum_pos_1_S)
law_poly_product_S_2 = law_convolution(law_sum_neg_1_S, law_sum_neg_2_S)
law_poly_product_S = law_convolution(law_poly_product_S_1, law_poly_product_S_2)
law_eT_S_prime = iter_law_convolution(law_poly_product_S, k)
print("  Distribution of e^T*S' calculated.")

# --- Calculate the distribution of s^T*E' ---
print("  --- 2b: Calculating distribution of s^T*E' ---")
law_E_prime_neg = {-val: prob for val, prob in law_E_prime.items()}
law_E_prime_2x = {2 * val: prob for val, prob in law_E_prime.items()}
law_E_prime_neg_2x = {-2 * val: prob for val, prob in law_E_prime.items()}
law_sum_pos_2_E = iter_law_convolution(law_E_prime_2x, num_twos)
law_sum_pos_1_E = iter_law_convolution(law_E_prime, num_ones)
law_sum_neg_1_E = iter_law_convolution(law_E_prime_neg, num_minus_ones)
law_sum_neg_2_E = iter_law_convolution(law_E_prime_neg_2x, num_minus_twos)
law_poly_product_E_1 = law_convolution(law_sum_pos_2_E, law_sum_pos_1_E)
law_poly_product_E_2 = law_convolution(law_sum_neg_1_E, law_sum_neg_2_E)
law_poly_product_E = law_convolution(law_poly_product_E_1, law_poly_product_E_2)
law_sT_E_prime = iter_law_convolution(law_poly_product_E, k)
print("  Distribution of s^T*E' calculated.\n")

# 3. Combine all distributions to get the final distribution of E_total
print("Step 3: Combine all distributions to get the final distribution of E_total...")
law_sT_E_prime_neg = {-val: prob for val, prob in law_sT_E_prime.items()}
law_intermediate = law_convolution(law_eT_S_prime, law_sT_E_prime_neg)
law_E_total = law_convolution(law_intermediate, law_e_double_prime)
print("Final noise distribution law_E_total calculated.\n")

# ==============================================================================
# Calculate the final Decryption Failure Rate (DFR)
# ==============================================================================
print("--- Result Analysis ---")
failure_threshold = (t_repetition * q) / 4.0
print(f"Decryption failure threshold: {failure_threshold}")
single_coeff_failure_prob = tail_probability(law_E_total, failure_threshold)
print(f"Failure probability for a single coefficient: {single_coeff_failure_prob}")
if single_coeff_failure_prob > 0:
    print(f"≈ 2^({math.log2(single_coeff_failure_prob):.2f})")

total_dfr = n * single_coeff_failure_prob
print(f"\nTotal Decryption Failure Rate (DFR) for the scheme: {total_dfr}")
if total_dfr > 0:
    print(f"≈ 2^({math.log2(total_dfr):.2f})")