# Protein-Shell Immobilized Enzymes in a Batch Reactor

A computational model of enzymes encapsulated in a porous protein shell and reused over many batch cycles. The model couples diffusion through the shell pores, Michaelis–Menten kinetics with product inhibition in the core, and enzyme deactivation, so we can ask: **which shell design keeps a pharmaceutical enzyme productive for the most cycles while still giving good conversion, yield, and product quality?**

*Project proposed by Dileep Reddy.*

## Team

| Name | GitHub | Role / focus |
|------|--------|--------------|
| Dileep Reddy | `@TODO` | Proposer; kinetics and deactivation |
| Alex | `@TODO` | TODO |
| TODO | `@TODO` | TODO |
| rwest | `@rwest` | Instructor |

## Plan

1. **Week 1 – Minimal model.** Single batch, one particle size, well-mixed core, linear shell resistance, Michaelis–Menten with competitive product inhibition, first-order deactivation. Reproduce the free-enzyme limit (no shell) as a sanity check.
2. **Week 2 – Repeated cycles.** Carry enzyme activity from cycle to cycle, add per-cycle shell degradation/leaching, compute conversion per cycle and cycles-to-50%-activity.
3. **Week 3 – Design study.** Sweep shell pore radius, shell thickness, and substrate size; map the trade-off between first-cycle conversion (thin, open shells) and lifetime (thick, protective shells).
4. **Week 4 – Extensions and write-up.** Temperature dependence (Arrhenius for both k<sub>cat</sub> and k<sub>d</sub>), mixing rate via an external film resistance, a side reaction for product purity, and final figures.

Each step should be a working, plotted result before moving to the next.

## Inputs, Physics, Outputs

| Inputs | Physics | Outputs |
|--------|---------|---------|
| Temperature and mixing rate | Reaction kinetics: Michaelis–Menten reaction in the core, with product inhibition | Conversion |
| Batch time and number of cycles | Mass transfer: diffusion of substrate into the shell and of product out | Enzyme activity per cycle |
| Shell pore radius and thickness | Enzyme deactivation | Effectiveness factor |
| Shell material | Shell degradation | Cycles to 50% activity |
| Substrate size and concentration | | Product yield and quality |

## Minimal Model

### Geometry
Spherical particles: an enzyme-loaded core of radius *R*<sub>c</sub> surrounded by a porous protein shell of thickness δ (outer radius *R*<sub>o</sub> = *R*<sub>c</sub> + δ). *N* particles are suspended in a liquid volume *V*<sub>L</sub>.

### Assumptions (minimal version)
- Bulk liquid and particle core are each well mixed (no radial profile inside the core).
- The shell is a pure diffusion barrier with no reaction and no accumulation (pseudo-steady flux).
- Isothermal; external film resistance neglected at first (added later to capture mixing rate).
- Only substrate S → product P; no side reactions yet.

### Hindered diffusion through the shell
Effective diffusivity depends on the ratio of solute radius to pore radius, λ = *r*<sub>s</sub>/*r*<sub>p</sub>, using the Renkin equation:

$$
D_\text{eff} = D_0\,\frac{\varepsilon}{\tau}\,(1-\lambda)^2\left(1 - 2.104\lambda + 2.09\lambda^3 - 0.95\lambda^5\right), \qquad \lambda < 1
$$

Diffusive conductance of a spherical shell:

$$
G = \frac{4\pi D_\text{eff} R_c R_o}{R_o - R_c}
$$

(Computed separately for substrate and product, since they can have different sizes.)

### Kinetics in the core
Michaelis–Menten with competitive product inhibition:

$$
v = \frac{k_\text{cat} E\, S_c}{K_m\left(1 + P_c/K_i\right) + S_c}
$$

### Balances
With core volume *V*<sub>c</sub> = (4/3)π*R*<sub>c</sub><sup>3</sup>:

$$
\begin{aligned}
V_c \frac{dS_c}{dt} &= G_S (S_b - S_c) - V_c\, v \\
V_c \frac{dP_c}{dt} &= V_c\, v - G_P (P_c - P_b) \\
\frac{dS_b}{dt} &= -\frac{N}{V_L} G_S (S_b - S_c) \\
\frac{dP_b}{dt} &= \frac{N}{V_L} G_P (P_c - P_b) \\
\frac{dE}{dt} &= -k_d(T)\, E
\end{aligned}
$$

### Between cycles
At the end of each batch the particles are recovered and put into fresh substrate solution. Enzyme activity carries over, multiplied by (1 − *f*<sub>leach</sub>) to represent loss from shell degradation. A later version can let *r*<sub>p</sub> grow with cycle number instead.

### Derived outputs
- **Conversion:** *X* = 1 − *S*<sub>b</sub>(*t*<sub>batch</sub>)/*S*<sub>0</sub>
- **Yield:** *P*<sub>b</sub>(*t*<sub>batch</sub>)/*S*<sub>0</sub> (equal to *X* until a side reaction is added)
- **Effectiveness factor:** η = *v*(*S*<sub>c</sub>, *P*<sub>c</sub>) / *v*(*S*<sub>b</sub>, *P*<sub>b</sub>)
- **Relative activity:** *E*<sub>n</sub>/*E*<sub>0</sub> at the start of cycle *n*
- **Cycles to 50% activity:** first *n* with *E*<sub>n</sub>/*E*<sub>0</sub> ≤ 0.5

### Starting parameters
Order-of-magnitude placeholders, to be replaced with literature values once we choose a specific enzyme and substrate.

| Parameter | Symbol | Placeholder | Units |
|-----------|--------|-------------|-------|
| Free diffusivity (small molecule) | *D*<sub>0</sub> | 1 × 10⁻⁹ | m²/s |
| Shell porosity / tortuosity | ε/τ | 0.3 | – |
| Michaelis constant | *K*<sub>m</sub> | 1 | mM |
| Inhibition constant | *K*<sub>i</sub> | 5 | mM |
| Turnover number | *k*<sub>cat</sub> | 100 | s⁻¹ |
| Deactivation constant | *k*<sub>d</sub> | 1 × 10⁻⁵ | s⁻¹ |
| Core radius | *R*<sub>c</sub> | 50 | µm |
| Shell thickness | δ | 1–20 | µm |
| Leaching fraction per cycle | *f*<sub>leach</sub> | 0.01 | – |

## First Plot

**Conversion at the end of each batch vs. cycle number**, one curve per shell thickness (or pore radius), with free enzyme as a reference line.

- **x-axis:** cycle number (1 to ~50)
- **y-axis:** end-of-batch conversion *X* (0 to 1)
- **Curves:** δ = 1, 5, 10, 20 µm, plus free enzyme
- **Marker:** horizontal line at 50% of first-cycle performance

We expect thin shells to start high and fall fast, while thick shells start lower (diffusion-limited, low η) but last longer. Where the curves cross is the design trade-off the project is about.

A companion plot of η vs. Thiele-like modulus for the same designs will help explain *why* each curve starts where it does.

## Repository Structure (proposed)

```
.
├── README.md
├── environment.yml        # conda environment (numpy, scipy, matplotlib, jupyter)
├── src/
│   ├── model.py           # ODE right-hand side, parameter dataclass
│   ├── diffusion.py       # Renkin hindrance, shell conductance
│   └── cycles.py          # multi-cycle driver
├── notebooks/
│   └── 01_first_plot.ipynb
├── figures/
└── tests/
    └── test_limits.py     # e.g. no shell → free-enzyme result
```

## Getting Started

```bash
git clone https://github.com/TODO/TODO.git
cd TODO
conda env create -f environment.yml
conda activate protein-shell
jupyter lab
```

## Contributing

Everyone commits. Work on a branch, open a pull request, and ask a teammate to review before merging into `main`. Keep notebooks lightweight and move reusable code into `src/`.

## References

1. Thiele, E. W. (1939). Relation between catalytic activity and size of particle. *Industrial & Engineering Chemistry*, 31(7), 916–920.
2. Renkin, E. M. (1954). Filtration, diffusion, and molecular sieving through porous cellulose membranes. *Journal of General Physiology*, 38(2), 225–243.
3. Engasser, J.-M., & Horvath, C. (1973). Effect of internal diffusion in heterogeneous enzyme systems: evaluation of true kinetic parameters and substrate diffusivity. *Journal of Theoretical Biology*, 42(1), 137–155.
4. Sheldon, R. A., & van Pelt, S. (2013). Enzyme immobilisation in biocatalysis: why, what and how. *Chemical Society Reviews*, 42(15), 6223–6235.
5. Bailey, J. E., & Ollis, D. F. (1986). *Biochemical Engineering Fundamentals* (2nd ed.). McGraw-Hill.
6. Fogler, H. S. *Elements of Chemical Reaction Engineering*. Pearson.

TODO: add a reference specific to protein-shell (e.g. protein cage or protein-polymer capsule) encapsulation and one with measured reuse-cycle data for our chosen enzyme.
