import sys
# Add the directory containing the estimator library to Python's search path
sys.path.append('/lattice-estimator')

from estimator import *
from sage.all import sqrt, log, oo, RR, factorial, round
from estimator.nd import NoiseDistribution

# ==============================================================================
# Wrapper class to adapt probability dictionary for the estimator 
# ==============================================================================
class ReconciliationDistribution(NoiseDistribution):
    """A wrapper class to adapt the probability dictionary for the estimator."""

    def __init__(self, dist: dict, n: int = None):
        self.dist = dist
        
        # Ensure probabilities sum to 1.0
        total_prob = sum(self.dist.values())
        if abs(total_prob - 1.0) > 1e-9:
             self.dist = {k: v/total_prob for k,v in self.dist.items()}
             
        mean = sum(val * prob for val, prob in self.dist.items())
        variance = sum(((val - mean)**2) * prob for val, prob in self.dist.items())

        super().__init__(
            n=n, 
            mean=RR(mean),
            stddev=RR(sqrt(variance)),
            bounds=(min(self.dist.keys()), max(self.dist.keys()))
        )

    def resize(self, new_n):
        return ReconciliationDistribution(dist=self.dist, n=new_n)

    def support_size(self, fraction=1.0):
        """
        Calculates the size of the search space (Multinomial coefficient) for fixed-weight distributions.
        Size = n! / (count_1! * count_2! * ... * count_m!)
        """
        # Numerator: n!
        total_perms = factorial(self.n)
        denominator = 1
        
        # Denominator: product of factorials of coefficient counts
        for prob in self.dist.values():
            count = int(round(prob * self.n))
            denominator *= factorial(count)
            
        total_size = total_perms / denominator
        
        return total_size * fraction

    def __call__(self, *args, **kwargs):
        return self.dist

    def __repr__(self):
        return f"ReconDist(μ={float(self.mean):.2f}, σ={float(self.stddev):.2f})"

# ==============================================================================
# Helper function: Build the reconciliation error distribution
# ==============================================================================
def build_reconciliation_law(t):
    # Step size is fixed to 2
    STEP = 2 
    
    # 1. Generate anchor points
    anchors = []
    val = 1
    while val < t:
        anchors.append(val)
        val += STEP
        
    counts = {}
    
    # 2. Iterate through the residue system of t
    for x in range(t):
        candidates = []
        for r in anchors:
            # Compute distance on the ring
            diff = (x - r + t//2) % t - t//2
            
            # Handle t/2 cases by generating both positive and negative directions
            if t % 2 == 0 and abs(diff) * 2 == t:
                candidates.append(abs(diff))  
                candidates.append(-abs(diff)) 
            else:
                candidates.append(diff)
        
        # Sort deterministically
        candidates.sort(key=lambda d: (abs(d), -d))
        
        # Select the best candidate deterministically based on sorting rules
        best = candidates[0]
        counts[best] = counts.get(best, 0.0) + 1.0

    # 3. Normalize
    dist_dict = {val: count/t for val, count in counts.items()}
    
    print(f"Debug [t={t}]: Error Dist = {dict(sorted(dist_dict.items()))}")
    return ReconciliationDistribution(dist=dist_dict)

print("="*60)
print(" Security Assessment for Lore PKE Scheme")
print("="*60)

# 2. Define the parameter sets for the Lore scheme based on Table 2 specifications
lore_params = {
    "Lore-128": {
        "n_ring": 512, "k": 1, "t": 2, 
        "s_dist_counts": {-2: 1, -1: 50, 0: 410, 1: 50, 2: 1}
    },
    "Lore-256": {
        "n_ring": 512, "k": 2, "t": 2, 
        "s_dist_counts": {-2: 12, -1: 140, 0: 720, 1: 140, 2: 12}
    },
    "Lore-384": {
        "n_ring": 512, "k": 3, "t": 4, 
        "s_dist_counts": {-2: 60, -1: 270, 0: 876, 1: 270, 2: 60}
    },
    "Lore-512": {
        "n_ring": 768, "k": 3, "t": 4, 
        "s_dist_counts": {-2: 60, -1: 420, 0: 1344, 1: 420, 2: 60}
    }
}

# 3. Loop through each parameter set and perform the evaluation
for name, params in lore_params.items():
    print(f"\n--- Evaluating parameter set: {name} ---\n")
    print(">>> (A) Key Recovery Attack (MLWE -> LWE)")

    t = params["t"]
    q = t * 257 
    RECONCILIATION_T = t 

    lwe_n = params["n_ring"] * params["k"]
    lwe_m = params["n_ring"] * params["k"]

    # Extract exact secret key distribution 
    s_dist_counts = params["s_dist_counts"]
    
    # Validation check for secret coefficient counts
    current_sum = sum(s_dist_counts.values())
    target_sum = lwe_n  # total coefficients = n * k
    
    if current_sum != target_sum:
        raise ValueError(
            f"Configuration Error: Sum of secret counts ({current_sum}) "
            f"does not match n_ring * k ({target_sum}) for parameter set '{name}'."
        )
        
    s_dist_prob_dict = {
        val: count / target_sum 
        for val, count in s_dist_counts.items()
    }
    dist_s = ReconciliationDistribution(dist=s_dist_prob_dict)

    # Generate error distribution
    if RECONCILIATION_T >= 2:
        dist_e = build_reconciliation_law(RECONCILIATION_T)
    else:
        dist_e = ND.Uniform(0, 1) 

    # Define the LWE problem instance
    mlwe_params = LWE.Parameters(
        n=lwe_n,
        q=q, 
        Xs=dist_s,
        Xe=dist_e,
        m=lwe_m,
        tag=f"{name}-KeyRecovery"
    )

    print("Secret distribution: Fixed-weight")
    print(f"Calculated -> {dist_s}")
    print(f"Error distribution: {dist_e}")
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