# Verification scope

Verified 25 September 2026 with existing shared environments; no package installation or upgrade was required.

- Ten standard-library tests passed: endmember stoichiometry, all recipe element balances, analytical fayalite yield, no hydrogen from the pure Mg hydration path, water limitation, activation-volume signs, effective stress, zero conversion, Gibbs pressure units and radiolysis energy units.
- All six reference reaction-path cases and the 162-row temperature/pressure/activation-volume sweep completed.
- Six PHREEQC cases each at 25, 80 and 95 °C completed with nonnegative solid amounts and cation relative-balance errors below 2e-14. The unsupported 100 °C input was rejected. These checks verify accounting and operation, not natural assemblages or reaction rates.
- The PHREEQC runner checks the normalized database hash; formula units for enstatite and ferrosilite are explicitly matched to the database.
- Reference numbers are retained under `reference/`. The equilibrium reference uses 0.1 kg rock, 1 kg initially dilute NaCl water and 1 atm. Reaction-path references use the separate conditions in `inputs/scenarios.json`.

The default calculations were tested with Python 3.11.16. Optional plotting was tested with NumPy 2.5.2 and Matplotlib 3.11.1 in an existing numerical environment. Optional equilibrium was tested with phreeqpy 0.6.0 and IPhreeqc 3.7.3 on macOS arm64. Package pins describe the tested optional Python layers; the native PHREEQC library and thermodynamic database must also be supplied.

No field calibration, full solid-solution equilibrium, coupled fracture propagation, or pressure-dependent mineral equation of state is claimed.
