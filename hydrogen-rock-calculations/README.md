# Hydrogen source rock calculations

Small, standalone calculation examples accompanying Section 1.5 of *Geological and economic potential of white hydrogen as a renewable resource* by Sander Bert da Costa. This folder can be used as a GitHub repository independently of the paper and the reservoir models. It has not been published to GitHub.

The calculations cover six declared ideal rock recipes: dunite, harzburgite, pyroxenite, mafic rock, magnetite-rich banded iron formation and a granite radiolysis host. These are illustrative inputs, not measured rock compositions or calibrated production forecasts.

## Run the reaction-path calculation

Use an existing Python 3.11 or newer interpreter. The default calculation and tests need only the standard library. Set `ROCK_PYTHON` to its absolute path; no environment activation is needed.

```sh
export ROCK_PYTHON="$HOME/.local/share/python-envs/darts-py311/bin/python"
"$ROCK_PYTHON" -B run_calculations.py
"$ROCK_PYTHON" -B -m unittest discover -s tests -v
```

For another computer, replace the interpreter path with a compatible local interpreter. Run these commands from this folder. Reuse a suitable shared environment rather than creating a new environment for this folder.

`inputs/scenarios.json` holds the mineral mass fractions and all scenario assumptions. The script writes per-mineral amounts, hydrogen and water balances, expansion and cracking indices, and a temperature/pressure/activation-volume sensitivity table into `results/`. Use `--input PATH` and `--output PATH` to select alternatives. Relative supplied paths are relative to the working directory; default paths are relative to the script.

With Matplotlib available, add `--plot` to produce the temperature/pressure figure. The optional plotting versions are recorded in `requirements-plot.txt`. Check shared-environment compatibility before installing or changing packages; use the selected absolute interpreter with `-m pip` if installation is necessary.

## What the calculation means

The reaction-path calculation applies balanced, prescribed endmember reactions to a finite initial inventory. An illustrative Arrhenius rate, accessible mineral fraction and water limit determine conversion. Mineral amounts give the change in solid volume. A constrained-expansion stress estimate and assumed flaw geometry give a fracture screening index. The index does not solve crack propagation or predict permeability.

Temperature and pressure affect the declared kinetics; they do not select stable mineral phases in this calculation. The Gibbs-energy helper accepts thermodynamic coefficients supplied by the caller and is not automatically used to infer mineral stability. Mineral densities are nominal reference values. Radiolysis uses deposited power and an illustrative yield separately from the water–rock reactions; it does not infer radiation power from granite composition.

The default dunite example gives **0.08244 g H₂ per kg initial rock**, **6.43% solid-volume increase**, and a **1.04 fracture index** after the specified 1000 years at 80 °C and 20 MPa. These values follow the input assumptions. They are not a comparison of natural rock productivity or a measured cracking threshold. See [the equations and worked examples](docs/calculations.md).

## Optional PHREEQC equilibrium example

`run_equilibrium.py` calculates a separate, restricted pure-phase equilibrium assemblage using `phreeqpy` and a working IPhreeqc native library. It is not required for the reaction-path calculation. The verified environment used Python 3.11.16, phreeqpy 0.6.0 and IPhreeqc 3.7.3 on macOS arm64. Other platforms need their own compatible IPhreeqc library; pass its absolute path with `--library` if automatic discovery is unsuitable.

Obtain the external `supcrtbl.dat` database from the pinned upstream location and check the provenance in [inputs/equilibrium_database.json](inputs/equilibrium_database.json). Keep the database outside this repository; it is not redistributed here. The runner verifies its normalized-text SHA-256 and rejects an untested replacement. Do not merely change the expected hash to bypass the check: phase formula units and thermodynamic conventions must first be validated.

```sh
"$ROCK_PYTHON" -B run_equilibrium.py --database /absolute/path/to/supcrtbl.dat
```

This example is restricted to **25–95 °C, 1 atm**, with 0.1 kg rock and 1 kg initially dilute NaCl water. The database tabulates constants along water's saturation-pressure curve, so this is a low-pressure approximation. Changing solution pressure alone does not provide consistent high-pressure mineral thermodynamics. It is not a brine or subsurface assemblage prediction.

Only the listed pure phases can form. Fe–Mg solid solutions and ferric serpentine are absent. H₂ starts with zero gas and can appear at unit fugacity; no hydrogen supply is imposed. This is a fixed-fugacity gas endmember, not a finite-headspace or gas-mixture pressure calculation. No kinetic time is assigned. The initial water/rock ratio and pressure differ from the reaction-path example, so their yields are not directly comparable.

Outputs include reusable PHREEQC input files, final phase amounts, aqueous and gas hydrogen, pH, database identity and cation-balance checks. `reference/` contains compact reference results for checking reproducibility. `results/` and native libraries are excluded from version control.

## Files and checks

- `rock_model/model.py`: stoichiometry, rate, Gibbs-energy helper, volume, fracture screening and radiolysis.
- `run_calculations.py`: recipe calculations and sensitivity tables.
- `run_equilibrium.py`: optional low-pressure equilibrium example.
- `inputs/`: explicit recipes and database provenance.
- `docs/`: equations, literature and verification scope.
- `tests/`: conservation, limiting cases and unit/sign checks.
- `reference/`: small results from the verified examples.

No DARTS installation is needed by these scripts. A reservoir model would still need to convert reaction extents to a generation rate per bulk volume and couple transport, chemistry and mechanics. This folder does not claim that coupling.

The code has not been assigned a redistribution licence by its author; see `NOTICE.md`. No external database or compiled runtime is bundled.
