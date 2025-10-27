# File: my-tests/test_lore_pke_params.py
# This script provides a security analysis for the Lore PKE scheme,
# with a correctly defined sparse secret distribution and a selectable error distribution.

# --- Code required for the solution ---
import sys
# Add the directory containing the estimator library to Python's search path
sys.path.append('/lattice-estimator')
# --- Code required for the solution ---

# 1. Import all functionalities from the estimator library
from estimator import *
# === Key Fix 1: Add RR to the import list ===
from sage.all import sqrt, log, oo, RR
# Import the base class
from estimator.nd import NoiseDistribution

# ==============================================================================
# Wrapper class to adapt our probability dictionary for the estimator
# ==============================================================================
class ReconciliationDistribution(NoiseDistribution):
    """A wrapper class to adapt our probability dictionary for the estimator."""

    # 1. Fix: Make __init__ accept n=None and pass it to the base class
    def __init__(self, dist: dict, n: int = None):
        """
        Manually defined constructor.
        :param dist: Stores our probability dictionary, e.g., {-1: 0.5, 1: 0.5}
        :param n: The dimension of the distribution (required by the estimator)
        """
        # Store the core probability dictionary
        self.dist = dist
        
        # Calculate mean and standard deviation from the dictionary
        mean = sum(val * prob for val, prob in self.dist.items())
        variance = sum(((val - mean)**2) * prob for val, prob in self.dist.items())
        
        # Call the parent class's constructor to set all required attributes
        super().__init__(
            n=n,  # <--- Fix 1: Pass 'n' to the base class
            mean=RR(mean),
            stddev=RR(sqrt(variance)),
            bounds=(min(self.dist.keys()), max(self.dist.keys()))
        )

    def resize(self, new_n):
        # 2. Fix: Return a *new* instance instead of 'self'
        #    and pass the required 'dist' and 'new_n'
        return ReconciliationDistribution(dist=self.dist, n=new_n)

    def __call__(self, *args, **kwargs):
        # Allows some parts of the estimator to call it like a function
        return self.dist

    def __repr__(self):
        # Custom print output for easy debugging
        return f"ReconDist(μ={float(self.mean):.2f}, σ={float(self.stddev):.2f})"

# ==============================================================================
# Helper function: Build the reconciliation error distribution based on the t value
# ==============================================================================
def build_reconciliation_law(t):
    """
    Builds the corresponding error distribution based on the provided t value and reconciliation rule.
    """
    if t == 3:
        dist_dict = {-1: 1/3, 0: 1/3, 1: 1/3}
    elif t == 7:
        dist_dict = {-2: 1/14, -1: 2/7, 0: 2/7, 1: 2/7, 2: 1/14}
    elif t == 13:
        dist_dict = {-2: 1/26, -1: 4/13, 0: 4/13, 1: 4/13, 2: 1/26}
    else:
        raise ValueError("Unsupported t value for reconciliation. Please choose 3, 7, or 13.")
    
    # Return an instance of our new wrapper class
    return ReconciliationDistribution(dist=dist_dict)

print("="*60)
print(" Security Assessment for Lore PKE Scheme (eta=1, Sparse Secrets)")
print("="*60)

# 2. Define the parameter sets for the Lore scheme
lore_params = {
    "Lore-128": {"n_ring": 256, "k": 2, "q": 3*257, "eta": 1, "e_agg_count": 1, "hw": 80},
    "Lore-192": {"n_ring": 256, "k": 3, "q": 7*257, "eta": 1, "e_agg_count": 2, "hw": 120},
    "Lore-256": {"n_ring": 256, "k": 4, "q": 13*257, "eta": 1, "e_agg_count": 4, "hw": 120}
}

# Parameter for the selectable distribution
RECONCILIATION_T = 3

# 3. Loop through each parameter set and perform the evaluation
for name, params in lore_params.items():
    print(f"\n--- Evaluating parameter set: {name} (hw={params['hw']}) ---\n")

    # A. Security assessment for Key Recovery Attack (MLWE -> LWE)
    print(">>> (A) Key Recovery Attack (MLWE -> LWE)")

    lwe_n = params["n_ring"] * params["k"]
    lwe_m = params["n_ring"] * params["k"]

    # a) Distribution of the secret key s (remains unchanged)
    p = m = params['hw'] // 2
    dist_s = ND.SparseTernary(p, m, n=params["n_ring"])

    # b) Distribution of the noise e (now selectable)
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

    print("Secret distribution:", dist_s)
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

print("="*60)
print(" All parameter sets evaluated")
print("="*60)