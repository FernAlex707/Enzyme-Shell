# Detailed Roadmap: Protein-Shell Immobilized ADH

*The one-page project plan is in the [README](../README.md). This file keeps the full step-by-step roadmap.*

A computational modeling course project simulating **alcohol dehydrogenase (ADH)** encapsulated in a porous protein shell. The model couples enzyme kinetics, diffusion through the shell, external mass transfer, temperature effects, and repeated batch cycles with washing and pore blockage.

---

## Table of Contents

- [Research Questions](#research-questions)
- [Scope and Assumptions](#scope-and-assumptions)
- [Model Development](#model-development)
- [Numerical Methods](#numerical-methods)
- [Parameters and Validation](#parameters-and-validation)
- [Simulation Studies](#simulation-studies)
- [Report and Presentation](#report-and-presentation)
- [Timeline](#timeline)
- [Practical Advice](#practical-advice)
- [Requirements](#requirements)
- [Suggested Repository Layout](#suggested-repository-layout)

---

## Research Questions

1. How do **pore size** and **pore number** change the effectiveness factor?
2. How fast does **activity decay** over repeated cycles?
3. What is the **optimal pore design**?

## Scope and Assumptions

**Phase 1 (Week 1)**

| Item | Choice |
|---|---|
| Reaction system | One system, e.g. ethanol + NAD⁺ ⇌ acetaldehyde + NADH, or acetaldehyde reduction |
| Reactor type | Batch with repeated cycles (supports washing and cycle-number studies) |
| Geometry | Spherical shell of radius *R*, enzyme inside, substrate diffusing in through pores |
| Spatial model | 1D radial |

**Initial assumptions**

- Isothermal (temperature added later)
- Pressure drop included
- One shell size per run
- Uniform enzyme concentration inside the shell

---

## Model Development

Build from simple to complex, testing each layer before adding the next.

### Step 1: Free-enzyme baseline

Michaelis–Menten kinetics:

$$v = \frac{V_{max} \, S}{K_m + S}$$

Batch ODEs:

$$\frac{dS}{dt} = -v, \qquad \frac{dP}{dt} = +v$$

Validate against the analytical (integrated) Michaelis–Menten solution.

### Step 2: Diffusion through the shell

Reaction–diffusion in spherical coordinates:

$$\frac{\partial C}{\partial t} = D_{eff} \, \frac{1}{r^2} \frac{\partial}{\partial r}\left(r^2 \frac{\partial C}{\partial r}\right) - v(C)$$

Effective diffusivity:

$$D_{eff} = \frac{D \, \varepsilon \, H}{\tau}$$

- **Porosity** ε is derived from pore number and pore size.
- **Hindrance factor** *H* follows the Renkin equation, based on the ratio of solute size to pore size.
- **Tortuosity** τ is a model parameter.

Compute the **Thiele modulus** (φ) and **effectiveness factor** (η). These are the key metrics.

### Step 3: External mass transfer (mixing rate)

Surface boundary condition:

$$\text{flux} = k_L \left(C_{bulk} - C_{surface}\right)$$

Link *k*<sub>L</sub> to mixing rate through a Sherwood correlation, Sh = *f*(Re, Sc).

### Step 4: Byproducts, intermediates, and cofactor

- ADH requires NAD(H): treat it as fixed or as a second species.
- Add product inhibition, and a competitive or ordered bi-bi kinetic form if feasible.
- Track each species with its own diffusivity.

### Step 5: Temperature

- Arrhenius dependence for *k*<sub>cat</sub>
- Thermal deactivation:

$$\frac{dE}{dt} = -k_d(T) \, E$$

### Step 6: Cycles, washing, and pore blockage

- **Washing:** between cycles, remove the bulk solution and optionally allow some enzyme to leach out (depending on pore size relative to enzyme size).
- **Pore blockage:** shrink the effective pore radius or open fraction with accumulated product or aggregate, e.g.

$$\frac{d\varepsilon}{dt} = -k_b \cdot (\text{foulant})$$

  This lowers *D*<sub>eff</sub> cycle by cycle.
- **Activity per cycle** = (intrinsic deactivation) × (diffusional limitation from blockage)

---

## Numerical Methods

**Phase 2 (Week 2)**

**Method:** Method of lines. Finite-volume discretization in *r* (about 50–200 cells), then integrate with `scipy.integrate.solve_ivp` using `BDF` or `Radau`, since the system is stiff.

**Verification checks:**

- Mesh and time-step convergence
- Mass balance closure
- Limit cases: φ → 0 gives η → 1; large φ gives η ≈ 1/φ

**Tools:** Python (NumPy, SciPy, Matplotlib). MATLAB is acceptable if the course uses it.

---

## Parameters and Validation


- Collect literature values: ADH *K*<sub>m</sub>, *k*<sub>cat</sub>, *V*<sub>max</sub>, diffusivities of substrate and NAD(H), and typical pore sizes for protein shells and cages.
- Fit unknown parameters (τ, *k*<sub>d</sub>, *k*<sub>b</sub>) with `scipy.optimize.least_squares`.
- Hold out some data for validation, or compare trends qualitatively if quantitative data is sparse.

---

## Simulation Studies

**Phase 3 (Week 3)**

| Vary | Observe |
|---|---|
| Pore size, pore number | Effectiveness factor, conversion, yield |
| Pore size distribution (lognormal vs. uniform) | Average *D*<sub>eff</sub>, blockage sensitivity |
| Shell size | Diffusion path length, η |
| Mixing rate | External resistance, activity |
| Temperature | Activity vs. stability trade-off |
| Cycle number | Activity decay curve, pore blockage |

- Run a **global sensitivity analysis** (Sobol or one-at-a-time) to rank which inputs matter most.
- Finish with a **design map**, e.g. a heatmap of conversion over pore size × pore number.

---

## Report and Presentation


**Structure:** introduction and motivation, governing equations, numerical method, verification, results, discussion (design recommendations, limitations), conclusion.

**Key figures:**

- Radial concentration profiles
- η vs. φ
- Activity vs. cycle number
- Sensitivity ranking
- Design heatmap

**Limitations to state clearly:** 1D geometry, lumped blockage model, literature-derived parameters.

---

## Timeline

| Week | Milestone |
|---|---|
|   | Scope, geometry, research questions, assumptions |
|  1 | Free-enzyme Michaelis–Menten model validated |
|   | Shell diffusion model running, η computed |
|   | External mass transfer, cofactor and byproducts added |
| 2 | Temperature, cycles, washing and pore blockage; numerics verified |
|   | Parameters collected and fitted |
|   | Parameter sweeps and sensitivity analysis |
| 3 | Design map, figures, draft report |
|   | Final report and presentation |

---

## Practical Advice

- Have a **working minimal model by Week 3**. Add physics one layer at a time and test each addition.
- Treat **pore size distribution as a stretch goal**. Do uniform pores first.
- Keep a **parameter table** with a source and units for every number.
- **Nondimensionalize early** (Thiele modulus, Biot number, Damköhler number).

---

## Requirements

- Python 3
- NumPy
- SciPy
- Matplotlib

```bash
pip install numpy scipy matplotlib
```

## Suggested Repository Layout

```
.
├── README.md
├── params/          # Parameter table (values, units, sources)
├── src/             # Model code (kinetics, diffusion, BCs, cycles)
├── tests/           # Verification: analytical MM, limit cases, mass balance
├── notebooks/       # Exploratory runs and parameter sweeps
├── figures/         # Generated plots
└── report/          # Draft and final report, presentation
```
