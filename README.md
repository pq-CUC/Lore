# Lore: An LWE-based Key Encapsulation Mechanism with Variable Modulus and CRT Compression

This repository contains the complete implementation and analysis tools for the **Lore** cryptosystem.

The scheme introduces two primary innovations to effectively balance security, communication bandwidth, and the Decryption Failure Rate (DFR):

1.  **Variable Modulus**: The scheme utilizes a composite modulus $t \cdot q$. The parameter $t$ is varied across different security levels to achieve a DFR consistent with the targeted security strength, while $q$ is fixed to allow for efficient Number Theoretic Transform (NTT) operations.
2.  **CRT Compression**: It leverages the Chinese Remainder Theorem (CRT) to compress public keys and ciphertexts. This technique, inspired by Learning with Rounding (LWR), reduces bandwidth by only storing part of the CRT representation and treating the resulting difference as cryptographic noise.

Secret keys in Lore are sampled from a **fixed-weight distribution** to further minimize the DFR and improve computational performance.

This repository includes:
* A C reference implementation of the Lore scheme.
* Python-based tools for security estimation against known lattice attacks.
* Python scripts for analyzing the Decryption Failure Rate (DFR).

## Directory Structure

```
.
├── Lore/
│   ├── Reference implementation/   # C reference implementation
│   ├── Security estimation/        # Python security and DFR analysis tools
│   └── Lwe_SummaCrypt-3.pdf        # The principle document for the Lore scheme
└── README.md
```

## 1. C Reference Implementation

Located in `Lore/Reference implementation/ref/`, this is a full C implementation of the Lore scheme, supporting three security levels (128-bit, 192-bit, and 256-bit) as defined in the paper.

### Compilation

The project uses a `Makefile` that defines different compilation flags for each security level (`-DLORE_LEVEL=1`, `-DLORE_LEVEL=2`, `-DLORE_LEVEL=3`).

**Compiling All Levels**

To compile all test programs for all security levels, run:

```bash
cd "Lore/Reference implementation/ref"
make
```

This will generate the following executables:
* `test_lore_1`, `test_lore_2`, `test_lore_3`
* `test_speed_1`, `test_speed_2`, `test_speed_3`
* `PQCgenKAT_pke_1`, `PQCgenKAT_pke_2`, `PQCgenKAT_pke_3`

**Compiling a Specific Level**

The `Makefile` is configured with separate targets for each program and level. To compile a specific program, simply specify its name. For example:

To compile only the functionality test for Level 1 (128-bit security):
```bash
make test_lore_1
```

To compile only the speed test for Level 2 (192-bit security):
```bash
make test_speed_2
```

### Usage

* **Functionality Test**: Run `./test_lore_1` to verify the correctness of key generation, encryption, and decryption.
* **Performance Benchmark**: Run `./test_speed_1` to measure the cycle counts for the scheme's core cryptographic operations.
* **Known Answer Tests (KAT)**: Run `./PQCgenKAT_pke_1` to generate KAT files, which are useful for ensuring compliance and correctness.

To clean up the build artifacts, run:
```bash
make clean
```

## 2. Security & DFR Analysis

The `Lore/Security estimation/` directory contains Python scripts for a thorough analysis of the scheme's parameters.

### Security Estimation

The `Lore_security/` subdirectory contains a lattice security estimator based on SageMath. It can be used to evaluate the scheme's hardness against the best-known primal and dual lattice attacks.

**Example Usage:**
These scripts assess the security of the defined Lore parameter sets. You can run them within a SageMath environment:

```bash
# Ensure you are in the Lore/Security estimation/Lore_security/ directory

# For 128-bit security level parameters:
sage ./my-lwe-tests/test_my_params_crt1.py

# For 192-bit and 256-bit security level parameters:
sage ./my-lwe-tests/test_my_params_crt2-r.py
```

The scripts will output the estimated classical and quantum security levels in bits for each parameter set.

### DFR Calculation

The `Lore_failure/` subdirectory contains a script to compute the scheme's Decryption Failure Rate.

**Example Usage (`Lore_failure_crt_e.py`):**
This script calculates the probability of a single coefficient failing decryption and aggregates it to determine the total DFR of the scheme.

```bash
# Ensure you are in the Lore/Security estimation/ directory
python ./Lore_failure/Lore_failure_crt_e.py
```

## Citation

If you use Lore or the tools in this repository in your work, please cite the accompanying paper:

> [Author(s)], "**Lore: An LWE-based Key Encapsulation Mechanism with Variable Modulus and CRT Compression**", [Conference/Journal, Year].