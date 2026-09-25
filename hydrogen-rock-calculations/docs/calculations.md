## 1.5 Calculations for source rocks and mineral assemblages

Hydrogen generation can be calculated at different levels. A mineral inventory gives the amount of iron available for oxidation. A reaction path gives the mineral amounts after a specified extent of alteration. An equilibrium calculation determines which of the allowed phases can coexist with the water. These calculations answer different questions. None gives a geological production rate without information on reaction speed, reactive surface and water access.

The calculation setup considers ideal dunite, harzburgite, pyroxenite, mafic rock and magnetite-rich banded iron formation (BIF). A granite example represents a possible radiolysis host. Their mineral proportions are declared inputs, rather than analyses of specific samples. Olivine and pyroxene are represented by separate Mg and Fe endmembers. This is useful for bookkeeping, but it does not reproduce the activities of natural solid solutions. The reaction-path results should therefore be read as conditional calculations. The separate PHREEQC example explores a restricted equilibrium assemblage.

### 1.5.1 Mineral amounts and hydrogen balance

For an initial rock mass $M_r$, mineral mass fraction $w_j$ and molar mass $M_j$, the amount of mineral $j$ is

$$n_{j,0}=M_r w_j/M_j. \tag{S1}$$

The prescribed endmember reactions are balanced for Mg, Fe, Si, O and H. Forsterite hydration changes solid volume without producing hydrogen:

$$2\,\mathrm{Mg_2SiO_4}+3\,\mathrm{H_2O}\rightarrow
\mathrm{Mg_3Si_2O_5(OH)_4}+\mathrm{Mg(OH)_2}. \tag{S2}$$

The Fe endmember can be represented by

$$3\,\mathrm{Fe_2SiO_4}+2\,\mathrm{H_2O}\rightarrow
2\,\mathrm{Fe_3O_4}+3\,\mathrm{SiO_2}+2\,\mathrm{H_2}. \tag{S3}$$

For the selected pyroxene paths,

$$3\,\mathrm{MgSiO_3}+2\,\mathrm{H_2O}\rightarrow
\mathrm{Mg_3Si_2O_5(OH)_4}+\mathrm{SiO_2}, \tag{S4}$$

$$3\,\mathrm{FeSiO_3}+\mathrm{H_2O}\rightarrow
\mathrm{Fe_3O_4}+3\,\mathrm{SiO_2}+\mathrm{H_2}. \tag{S5}$$

A possible magnetite oxidation path is

$$2\,\mathrm{Fe_3O_4}+\mathrm{H_2O}\rightarrow
3\,\mathrm{Fe_2O_3}+\mathrm{H_2}. \tag{S6}$$

Equation S6 is a limiting reaction, not evidence that all magnetite in BIF will react by this route. The code uses hematite as the ferric endpoint. Experiments by Geymond et al. (2023) identified partial conversion to metastable maghemite, so phase identity and reaction pathway must remain explicit. Oxidation by atmospheric oxygen does not have the same hydrogen balance as reduction of water. The mafic example holds anorthite and diopside inert in the simple path; it is not a complete basalt alteration model.

With reaction extent $\xi_r$ and signed stoichiometric coefficients $\nu_{jr}$, final amounts satisfy

$$n_j=n_{j,0}+\sum_r\nu_{jr}\xi_r,\qquad n_j\geq0. \tag{S7}$$

The code checks every elemental balance. When the proposed extents require more water than is present, all extents are reduced by the same limiting multiplier. Newly formed magnetite is not oxidised again in the olivine and pyroxene paths. This prevents a second hydrogen-producing reaction from being introduced without being declared. A separate general electron balance is $n_{H_2}=\Delta n_{Fe(III)}/2$ when oxidation of Fe(II) is the only electron donor and reduction of water is the only electron sink. Iron entering magnetite is only partly oxidised, whereas Fe can also remain ferrous in alteration minerals.

### 1.5.2 Temperature and pressure effects

The thermodynamic driving force is the reaction Gibbs energy at the actual activities. For dimensionless reaction quotient $Q$, gas activities use fugacity divided by its standard-state fugacity:

$$\Delta_rG(T,p)=\Delta_rG^\circ(T,p)+RT\ln Q,\qquad
K=\exp[-\Delta_rG^\circ/(RT)]. \tag{S8}$$

At equilibrium $Q=K$. A negative reaction Gibbs energy indicates a favourable forward direction for the specified reaction and state, but does not set its rate. Increasing hydrogen fugacity raises $Q$ for hydrogen-producing reactions and can reduce the driving force. Water activity and mineral activities also matter. This is why temperature, pressure and mineral proportions cannot be treated as independent universal production multipliers.

A local approximation, with constant reaction enthalpy, entropy and volume, is

$$\Delta_rG^\circ(T,p)\simeq\Delta_rH^\circ_{ref}
-T\Delta_rS^\circ_{ref}+\Delta_rV^\circ(p-p_{ref}). \tag{S9}$$

Here $\Delta_rS^\circ_{ref}=(\Delta_rH^\circ_{ref}-\Delta_rG^\circ_{ref})/T_{ref}$. Pressure is in Pa and reaction volume in m³/mol, giving J/mol. This approximation is implemented as a calculation function, but no unverified Gibbs or enthalpy values are assigned to the rock recipes. Larger temperature and pressure ranges need consistent heat-capacity, fluid and mineral-volume data. A one-bar stability temperature cannot be converted into a depth limit using the geothermal gradient alone.

The simple reaction-path calculation uses an illustrative first-order rate:

$$k(T,p)=k_{ref}\exp\!\left[-\frac{E_a}{R}\left(\frac1T-\frac1{T_{ref}}\right)
-\frac{\Delta V^\ddagger(p-p_{ref})}{RT}\right],\quad
k_{ref}=\frac{\ln2}{t_{1/2,ref}}. \tag{S10}$$

The activation volume $\Delta V^\ddagger$ controls the pressure dependence of the rate; it is not the reaction volume in Equation S9. Its sign is tested rather than assumed. For each accessible primary mineral, the reacted fraction before water limitation is $f_a[1-\exp(-kt)]$. The reference uses $T_{ref}=80$ °C, $p_{ref}=20$ MPa, $t_{1/2,ref}=1000$ years, $E_a=60$ kJ/mol and $f_a=0.25$. These values define a sensitivity example, not an experimental fit. The sweep uses 80, 150 and 250 °C; 10, 20 and 30 MPa; and activation volumes of −10, 0 and +10 cm³/mol. It does not select the stable minerals or impose a thermodynamic affinity term.

Published equilibrium studies show that Fe partitioning and secondary mineral stability change with temperature and initial composition (McCollom and Bach, 2009; Klein et al., 2013). The Arrhenius curves must therefore not be read as proof that the highest temperature always produces the most hydrogen. A measured rate law would additionally need reactive area, pH, activities and a driving-force term consistent with the selected reaction.

### 1.5.3 Expansion and a cracking calculation

Reaction-induced expansion is calculated from the initial and final mineral amounts. With nominal mineral density $\rho_j$,

$$V_s=\sum_j n_jM_j/\rho_j,\qquad
\epsilon_{v,rxn}=\frac{V_{s,final}-V_{s,initial}}{V_{b,0}},\qquad
V_{b,0}=\frac{V_{s,initial}}{1-\phi_0}. \tag{S11}$$

The densities are reference estimates for the endmembers, not pressure-dependent equations of state. For example, the forsterite and lizardite values are 3270 and 2550 kg/m³, close to the values reported in the Handbook of Mineralogy. Water incorporation increases solid mass; this must not be mistaken for a loss of mass conservation. Solid-volume change also differs from the total volume change of solids plus fluid.

For a simple mechanical screening, a fraction $\eta$ of initial pore space can accommodate the new solids. A constrained expansion stress scale is then

$$\sigma_{rxn}=K_{eff}\max(\epsilon_{v,rxn}-\eta\phi_0,0),\qquad
\sigma_3'=\sigma_3-\alpha_Bp_f. \tag{S12}$$

Compression is positive in the effective-stress expression. The trial tensile loading near an assumed flaw is $\sigma_{open}=\max(\sigma_{rxn}-\sigma_3',0)$. For flaw half-length $a$ and geometry factor $Y$,

$$K_I=Y\sigma_{open}\sqrt{\pi a},\qquad
I_{crack}=K_I/K_{IC}. \tag{S13}$$

An index of at least one is a positive screening result for the assumed flaw and loading, not a solved fracture prediction. The reference assumes $K_{eff}=1$ GPa, $\phi_0=0.05$, $\eta=0.25$, $\sigma_3=50$ MPa, $\alpha_B=1$, $a=1$ mm, $Y=1$ and $K_{IC}=1$ MPa√m. These are scenario inputs. Creep, dissolution, heterogeneous stress, crack interaction and fluid-pressure evolution are omitted. Calculating a permeability increase or crack aperture from this index would require an additional mechanical model.

Reaction-driven cracking may expose fresh mineral and admit water, but precipitation can also reduce permeability. Kelemen and Hirth (2012) and Evans et al. (2020) explain this competition in more complete mechanical models. Ross et al. (2025, corrected 2026) observed spatially variable alteration under flow; their results do not establish a universal permeability-increase law. The present screening keeps expansion, hydrogen generation and connected flow as separate outputs.

### 1.5.4 An equilibrium calculation for assemblages

The second script uses PHREEQC to equilibrate a finite mineral inventory with water. In thermodynamic terms, equilibrium minimises total Gibbs energy subject to conserved element amounts and nonnegative phase amounts:

$$\min_{\boldsymbol n\geq0}G(T,p,\boldsymbol n)
\quad\mathrm{subject\ to}\quad\mathbf A\boldsymbol n=\boldsymbol b. \tag{S14}$$

PHREEQC solves the corresponding mass-action and balance equations. The candidate list includes forsterite, fayalite, enstatite, ferrosilite, lizardite, brucite, magnetite, hematite, quartz, talc, greenalite, minnesotaite, diopside, anorthite, albite, microcline, annite and clinochlore. Phases initially absent can precipitate. Candidate phases excluded from the list cannot appear, even if they would be important in a natural rock. Pure Mg and Fe phases do not reproduce Fe–Mg solid-solution mixing or ferric serpentine chemistry.

The working example contains 0.1 kg rock and 1 kg water initially containing 0.001 mol/kg NaCl. It is restricted to 25–95 °C and 1 atm. Hydrogen gas starts at zero amount and can appear at unit fugacity; there is no imposed hydrogen supply. Final solid amounts, dissolved hydrogen, gas hydrogen, pH and cation balances are exported. This equilibrium endpoint has no elapsed reaction time and is not directly comparable to the rate calculation at 20 MPa.

The external SUPCRTBL-derived database records its equilibrium constants along the water saturation-pressure curve. Simply changing the PHREEQC solution pressure would not provide complete high-pressure mineral corrections. The script therefore checks the database identity and limits its own example to low pressure. Extension to subsurface pressure requires a consistent pressure-specific database and validation of activities, volumes and relevant solid solutions (Parkhurst and Appelo, 2013; Zhang et al., 2020).

### 1.5.5 Calculated examples and interpretation

For the reaction-path reference, all rocks start with 1 kg solid and 0.5 kg water. After one assumed half-time, 12.5% of each selected primary mineral has reacted, because only 25% was assigned as accessible. Water is sufficient in these examples. The mineral amounts are saved separately, so the hydrogen and expansion calculations can be traced back to the selected reactions.

| Ideal rock | H₂ g per kg rock | Solid volume change % | Cracking index |
|:--|--:|--:|--:|
| Dunite | 0.0824 | 6.43 | 1.04 |
| Harzburgite | 0.0768 | 5.84 | 0.73 |
| Pyroxenite | 0.0955 | 4.55 | 0.04 |
| Mafic rock | 0.0731 | 1.06 | 0.00 |
| Magnetite BIF | 0.2721 | 0.08 | 0.00 |
| Granite radiolysis host | 0.0000 | 0.00 | 0.00 |

The recipes intentionally differ in Fe content. The table is not a ranking of real rock productivity: it prescribes the same illustrative conversion law and different reaction endpoints. The granite row has no prescribed water–rock reaction; zero here does not mean that real granite cannot generate hydrogen. Radiolysis is calculated separately.

For the dunite example, 0.1 kg initial fayalite corresponds to 0.4907 mol. Conversion of 12.5% gives 0.04090 mol hydrogen, or 0.08244 g, from Equation S3. The Mg endmember produces no hydrogen but supplies most of the hydration-related expansion. The calculated solid-volume increase is 6.43%, and the bulk reaction strain is 0.0611. After subtracting the assigned pore allowance, the trial expansion stress is 48.6 MPa. Against 30 MPa effective confinement, the assumed flaw gives $K_I=1.04$ MPa√m. Its proximity to the chosen toughness shows why the result is sensitive to mechanical inputs rather than demonstrating that a real dunite will fracture.

The equilibrium example gives a different illustration. At 80 °C and 1 atm, the dunite recipe retains lizardite, brucite and magnetite among the selected solid phases and contains about 0.03272 mol total hydrogen for 0.1 kg initial rock. The pyroxenite recipe instead gives lizardite, talc and greenalite, with about $5.00\times10^{-8}$ mol hydrogen. Within this restricted model, Fe retained in secondary silicates changes the hydrogen outcome substantially. The amount depends on the permitted phases and the database; it is not a validated natural yield or a rate.

### 1.5.6 Radiolysis and connection to reservoir models

Radiolysis is driven by deposited radiation energy, rather than by the serpentinization activation energy. If $P_{abs}$ is power absorbed by water and $G_{H_2}$ is the yield in molecules per 100 eV,

$$n_{H_2,gross}=\frac{G_{H_2}P_{abs}t}{100\,e\,N_A}. \tag{S15}$$

Here $e=1.602176634\times10^{-19}$ J/eV, $t$ is in seconds and $N_A$ is Avogadro's constant. The example uses an illustrative low-LET yield of 0.45 molecules per 100 eV and deposited powers of $10^{-13}$ to $10^{-11}$ W per kg rock. Absorbed power is an input, not inferred from the mineral recipe. U, Th and K concentrations, decay energies, radiation type, geometry and energy transfer to water are needed to estimate it. Recombination and chemical consumption make gross generation different from accumulated hydrogen (Le Caër, 2011). Radiolysis does not, by itself, define a serpentine assemblage or hydration-expansion law.

The scripts are stored in the separate hydrogen-rock-calculations folder. They write mineral amounts, water consumption, hydrogen amounts, temperature-pressure sensitivities and mechanical screening values. Reaction balances, water limitation, limiting yields and pressure-unit conversions are checked automatically. The PHREEQC calculation additionally checks the cations remaining in water and solids.

A reservoir calculation would need generation per unit bulk volume and time, for example $q_{H_2}=M_{H_2}\dot n_{H_2}/V_b$, together with water consumption and retained solid mass. These quantities can inform a source term for DARTS, but the small scripts do not automatically replace the supply assumptions in Chapter 5. That step requires a rate law and geochemical domain appropriate to the reservoir. The calculations clarify the link between mineral composition, alteration and possible hydrogen supply without treating a maximum inventory or an equilibrium endpoint as a sustainable production rate.
