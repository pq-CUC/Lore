## Lattice Estimator (Third-Party Tool)

This directory contains a copy of the **Lattice Estimator**, an open-source tool for estimating the security of lattice-based cryptographic schemes. This tool is used by the parent project (**Lore**) for security analysis.

* **Original Repository**: `https://github.com/malb/lattice-estimator`
* **Original Authors & Contributors**: Please see the `README.rst` file in the `Lore_security` subdirectory.
* **License**: The Lattice Estimator is licensed under the **LGPLv3+ license**.

All credit for the **core Lattice Estimator tool** goes to its original authors. The estimator is included here for convenience to reproduce the security analysis results for the Lore scheme.

### Contributions to this Directory

The test scripts located at `Lore_security/my-lwe-tests/` were written to specifically evaluate the security parameters of the Lore scheme using this estimator.