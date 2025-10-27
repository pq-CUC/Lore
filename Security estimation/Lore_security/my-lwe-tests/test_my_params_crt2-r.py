# File: my-tests/test_lore_pke_params.py
# This script provides a security analysis for the Lore PKE scheme,
# using a Gaussian approximation for a user-definable, fixed-weight 5-value secret distribution.

# --- Code required for the solution ---
import sys

# Add the directory containing the estimator library to Python's search path
sys.path.append('/lattice-estimator')
# --- Code required for the solution ---

# 1. Import all functionalities from the estimator library
from estimator import *
# Import required SageMath functions and the base class
from sage.all import sqrt, log, oo, RR
from estimator.nd import NoiseDistribution


# ==============================================================================
# Wrapper class to adapt our probability dictionary for the estimator (for noise e)
# ==============================================================================
class ReconciliationDistribution(NoiseDistribution):
    """A wrapper class to adapt our probability dictionary for the estimator."""

    # 1. Fix: Make __init__ accept n=None
    def __init__(self, dist: dict, n: int = None):
        self.dist = dist
        mean = sum(val * prob for val, prob in self.dist.items())
        variance = sum(((val - mean) ** 2) * prob for val, prob in self.dist.items())

        super().__init__(
            n=n,  # 2. Fix: 'n' is now a known parameter
            mean=RR(mean),
            stddev=RR(sqrt(variance)),
            bounds=(min(self.dist.keys()), max(self.dist.keys()))
        )

    def resize(self, new_n):
        # 3. Fix: Return a *new* instance and pass new_n into it
        return ReconciliationDistribution(dist=self.dist, n=new_n)

    def __call__(self, *args, **kwargs):
        return self.dist

    def __repr__(self):
        return f"ReconDist(μ={float(self.mean):.2f}, σ={float(self.stddev):.2f})"


# ==============================================================================
# Helper function: Build the reconciliation error distribution based on the t value (for noise e)
# ==============================================================================
def build_reconciliation_law(t):
    if t == 3:
        dist_dict = {-1: 1 / 3, 0: 1 / 3, 1: 1 / 3}
    elif t == 7:
        dist_dict = {-2: 1 / 14, -1: 2 / 7, 0: 2 / 7, 1: 2 / 7, 2: 1 / 14}
    elif t == 13:
        dist_dict = {-2: 1 / 26, -1: 4 / 13, 0: 4 / 13, 1: 4 / 13, 2: 1 / 26}
    else:
        raise ValueError("Unsupported t value for reconciliation. Please choose 3, 7, or 13.")

    return ReconciliationDistribution(dist=dist_dict)


print("=" * 60)
print(" Security Assessment for Lore PKE Scheme (Custom 5-value Secret)")
print("=" * 60)

# 2. Define the parameter sets for the Lore scheme
lore_params = {
    "Lore-128": {"n_ring": 256, "k": 2, "q": 3 * 257, "eta": 1, "e_agg_count": 1},
    "Lore-192": {"n_ring": 256, "k": 3, "q": 7 * 257, "eta": 1, "e_agg_count": 1},
    "Lore-256": {"n_ring": 256, "k": 4, "q": 13 * 257, "eta": 1, "e_agg_count": 1}
}

# Parameter for the selectable distribution (for noise e)
RECONCILIATION_T = 13

# 3. Loop through each parameter set and perform the evaluation
for name, params in lore_params.items():
    print(f"\n--- Evaluating parameter set: {name} ---\n")

    # A. Security assessment for Key Recovery Attack (MLWE -> LWE)
    print(">>> (A) Key Recovery Attack (MLWE -> LWE)")

    lwe_n = params["n_ring"] * params["k"]
    lwe_m = params["n_ring"] * params["k"]

    # === Customize the count of each coefficient in the secret key s here ===
    # Note: The sum of all counts must be equal to n_ring (e.g., 256)
    s_dist_counts = {
        -2: 22,
        -1: 40,
        0: 132,
        1: 40,
        2: 22
    }
    # Verify if the total sum is correct
    if sum(s_dist_counts.values()) != params["n_ring"]:
        raise ValueError(
            f"Custom secret key coefficient sum {sum(s_dist_counts.values())} does not match n_ring {params['n_ring']}!")

    # a) Automatically calculate the standard deviation based on the custom parameters above and approximate with a Gaussian distribution
    s_mean = sum(val * count for val, count in s_dist_counts.items()) / params["n_ring"]
    s_variance = sum(((val - s_mean) ** 2) * count for val, count in s_dist_counts.items()) / params["n_ring"]
    secret_stddev = sqrt(s_variance)
    dist_s = ND.DiscreteGaussian(secret_stddev)

    # b) The distribution of noise e is determined by the value of RECONCILIATION_T (this logic remains unchanged)
    if RECONCILIATION_T in [3, 7, 13]:
        dist_e = build_reconciliation_law(RECONCILIATION_T)
    else:
        error_variance = params["e_agg_count"] * (params["eta"] / 2.0)
        error_stddev = sqrt(error_variance)
        dist_e = ND.DiscreteGaussian(error_stddev)

    # Define the LWE problem instance
    mlwe_params = LWE.Parameters(
        n=lwe_n,
        q=params["q"],
        Xs=dist_s,
        Xe=dist_e,
        m=lwe_m,
        tag=f"{name}-KeyRecovery"
    )

    print("Secret distribution: Gaussian Approx. of custom fixed-weight dist")
    # === Key fix: Corrected the f-string syntax ===
    print(
        f"       Custom counts per poly: -2: {s_dist_counts[-2]}, -1: {s_dist_counts[-1]}, 0: {s_dist_counts[0]}, 1: {s_dist_counts[1]}, 2: {s_dist_counts[2]}")
    print(f"       Calculated -> {dist_s}")
    if RECONCILIATION_T in [3, 7, 13]:
        print(f"Error distribution: Reconciliation law for t={RECONCILIATION_T} -> {dist_e}")
    else:
        print("Error distribution (Gaussian):", dist_e)
    print("LWE equivalent parameters:", mlwe_params)

    # --- Classical Computer Security Assessment ---
    print("\n--- [Classical] Key Recovery Evaluation Results ---")
    lwe_results_classical = LWE.estimate(mlwe_params, quiet=True, red_cost_model=RC.MATZOV)
    min_lwe_classical = oo
    for attack, cost in lwe_results_classical.items():
        print(f"Attack '{attack}': {cost!r}")
        if cost is not None and cost.get("rop") is not None:
            current_rop = cost.get("rop", oo)
            if current_rop != oo and current_rop < min_lwe_classical:
                min_lwe_classical = current_rop

    if min_lwe_classical != oo:
        print(f"==> [Classical] Key Recovery Security Level: {float(log(min_lwe_classical, 2)):.2f} bits")
    else:
        print("==> [Classical] Key Recovery Security Level: Not determined")

    # --- Quantum Computer Security Assessment ---
    print("\n--- [Quantum] Key Recovery Evaluation Results ---")
    lwe_results_quantum = LWE.estimate(mlwe_params, quiet=True, red_cost_model=RC.ChaLoy21)
    min_lwe_quantum = oo
    for attack, cost in lwe_results_quantum.items():
        print(f"Attack '{attack}': {cost!r}")
        if cost is not None and cost.get("rop") is not None:
            current_rop = cost.get("rop", oo)
            if current_rop != oo and current_rop < min_lwe_quantum:
                min_lwe_quantum = current_rop

    if min_lwe_quantum != oo:
        print(f"==> [Quantum] Key Recovery Security Level: {float(log(min_lwe_quantum, 2)):.2f} bits\n")
    else:
        print("==> [Quantum] Key Recovery Security Level: Not determined\n")

print("=" * 60)
print(" All parameter sets evaluated")
print("=" * 60)