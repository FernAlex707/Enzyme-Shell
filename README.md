# Computational Model of Protein-Shell Immobilized ADH

*One-page project plan. The full step-by-step roadmap and weekly timeline are in [docs/ROADMAP.md](docs/ROADMAP.md).*

## The system and the question

We model alcohol dehydrogenase (ADH) encapsulated in a porous protein shell used as a reusable catalyst in a repeated-batch reactor. Our reference case is AdhD from *Pyrococcus furiosus* packaged inside the bacteriophage P22 capsid, because it is the best-characterized ADH-in-a-shell system. We want to find out how pore size and pore number set the shell's diffusional resistance and effectiveness factor, how activity decays over repeated cycles through enzyme deactivation, leaching and pore blockage, and which pore design gives the best conversion over a run.

## Inputs → physics → outputs

| Inputs (what we vary) | Physics (what connects them) | Outputs (what we predict, and what we compare with) |
|---|---|---|
| Pore radius and pore number (→ porosity ε) | Michaelis–Menten kinetics, later with NAD(H) as a second substrate and product inhibition | Effectiveness factor η vs. Thiele modulus φ / Damköhler number Da |
| Shell radius *R* | Reaction–diffusion in a sphere: $\partial C/\partial t = D\,\frac{1}{r^2}\partial_r(r^2\partial_r C) - v(C)$ | Conversion and yield per batch |
| Enzyme loading per shell | Shell as a permeable boundary: $\text{flux} = P_{shell}(C_{out} - C_{in})$, with $P_{shell} = \varepsilon H D / (\tau L_{shell})$ and *H* from the Renkin equation | Radial concentration profiles |
| Mixing rate | External film: Sherwood correlation for $k_L$ (Sh → 2 for nanoparticles) | Activity vs. cycle number |
| Temperature | Arrhenius $k_{cat}(T)$; thermal deactivation $dE/dt = -k_d(T)E$ | Design map: conversion over pore size × pore number |
| Number of cycles | Washing step with size-dependent leaching; pore blockage $d\varepsilon/dt = -k_b\cdot(\text{foulant})$ | **Compare with:** free vs. encapsulated kinetics (Ref. 3), size-exclusion trends across P22 morphologies (Ref. 2), limit cases φ → 0 (η → 1) and large φ (η ≈ 1/φ) |

## Uncertain inputs

Rough first list. Milestone 2's sensitivity analysis starts here. All ranges still need checking against the original papers.

| Input | Range | Where the range came from |
|---|---|---|
| Pore radius | ~1.25–5 nm (2.5–10 nm diameter) | ~2.5 nm pores in P22 procapsid (Ref. 3); ~10 nm openings in the "wiffleball" form (Ref. 2) |
| Shell permeability $P_{shell}$ | 10⁻⁷ to 10 cm/s | Low end: values explored in the microcompartment model (Ref. 1); high end: Ref. 1's estimate for free diffusion through the open pore area. The gap is attributed to pore electrostatics, which Renkin ignores |
| *k*<sub>cat</sub> | 0.46–0.80 s⁻¹ | Encapsulated vs. free AdhD (Ref. 3, SI) |
| *K*<sub>m</sub> | 3.0–4.7 mM | Free vs. encapsulated AdhD (Ref. 3, SI) |
| Enzymes per shell | ~50–250 | ~204 AdhD per capsid (Ref. 2); fewer when co-assembled with scaffold protein (Ref. 3). Range is our guess |
| Shell radius *R* | 25–35 nm | 58 nm P22 diameter (Ref. 3); expanded forms are larger |
| Substrate *D* (small alcohol/ketone) | 0.8–1.3 × 10⁻⁹ m²/s | Typical small-solute values in water (ethanol ≈ 1.2 × 10⁻⁹); to confirm for our substrate |
| NAD(H) *D* | 2–7 × 10⁻¹⁰ m²/s | Stokes–Einstein estimate from molecular size; needs a literature value |
| Tortuosity τ | 1–3 | Assumed; a thin monolayer shell is likely close to 1. To be fitted |
| *k*<sub>d</sub> | 0 to ~9 × 10⁻⁵ s⁻¹ | Upper bound from AdhD's 130 min half-life at 100 °C (Machielsen et al. 2006); much lower at our temperatures |
| *k*<sub>b</sub> (pore blockage) | Unknown; span 0 to "50 % porosity loss in 10 cycles" | No source yet. To be fitted or treated as a scenario parameter |

## Three key references

1. **Jakobson, Tullman-Ercek, Slininger & Mangan (2017).** A systems-level model reveals that 1,2-propanediol utilization microcompartments enhance pathway flux through intermediate sequestration. *PLoS Comput Biol* 13: e1005525. [doi:10.1371/journal.pcbi.1005525](https://doi.org/10.1371/journal.pcbi.1005525). **Gives us the model:** a spherical reaction–diffusion model with Michaelis–Menten enzymes inside a protein shell and a permeability boundary condition, plus a parameter table, characteristic-time analysis and public MATLAB code.
2. **Selivanovitch, LaFrance & Douglas (2021).** Molecular exclusion limits for diffusion across a porous capsid. *Nat Commun* 12: 2903. [doi:10.1038/s41467-021-23200-1](https://doi.org/10.1038/s41467-021-23200-1). **Gives us pore data and a mechanism:** pore sizes for P22 morphologies, enzyme loading per capsid, measured size and charge cut-offs for substrates reaching encapsulated AdhD, and evidence that electrostatics and capsid "breathing" matter beyond a purely geometric pore model.
3. **Sharma & Douglas (2020).** Tuning the catalytic properties of P22 nanoreactors through compositional control. *Nanoscale* 12: 336–346. [doi:10.1039/C9NR08348K](https://doi.org/10.1039/C9NR08348K). **Gives us parameter values and data to compare with:** *k*<sub>cat</sub> and *K*<sub>m</sub> for free and encapsulated AdhD, and activity vs. packing density.

## The first thing we'll compute

Steady-state effectiveness factor for one substrate and one Michaelis–Menten enzyme inside a sphere of radius *R*, with the shell as a permeable boundary. First the well-mixed-interior version, which is a single algebraic equation:

$$\frac{3P_{shell}}{R}\,(C_{out} - C_{in}) = \frac{V_{max}\,C_{in}}{K_m + C_{in}}, \qquad \eta = \frac{v(C_{in})}{v(C_{out})}$$

Then the same problem with radial finite volumes, checked against it. We sweep $P_{shell}$ across the full uncertain range and plot η against $\text{Da} = V_{max}R / (3P_{shell}K_m)$. A back-of-envelope estimate with the values above puts Da ≈ 1 near $P_{shell}$ ≈ 5 × 10⁻⁷ cm/s, inside our uncertain range. Whether the shell limits the rate at all is therefore the first thing this model tells us. Diffusion *inside* the ~25 nm lumen looks negligible (internal φ² ~ 10⁻⁷), so the shell boundary is where the resistance sits.

## Who does what, for now

| Person | Responsibility |
|---|---|
| [Name] | Kinetics module (Michaelis–Menten → bi-bi, inhibition), free-enzyme validation |
| [Name] | Shell transport (Renkin, permeability, radial solver), numerical verification |
| [Name] | Literature and parameter table (values, units, sources), fitting |
| [Name] | Cycles, washing, blockage; sensitivity analysis and figures |
| Everyone | Weekly check of the parameter table and this README |

## Use of AI

We add a row every time we use an AI tool for this project.

| Date | Tool | What we used it for | How we checked it |
|---|---|---|---|
| Nov 2026 | Claude (Anthropic) | Converted our roadmap PDF into a README | Compared section by section against the PDF |
| Nov 2026 | Claude (Anthropic), with web search | Restructured the README into this plan; proposed the three references and first-pass parameter ranges; suggested treating the shell as a permeable boundary | **To do:** open each reference and confirm authors, title and every number we cite; redo the Da and φ estimates by hand |
