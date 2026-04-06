# ATLAS Di-Higgs Discrimination Project (Particle Physics)

**University of Manchester — MPhys Project, Sep 2019 – Jan 2020**  
**Supervisor:** Prof. Terry Wyatt  

---

**Full project report:** See [`report.pdf`](./report.pdf) for the complete methodology, results, and ROC curve analysis.

---

## Overview

This project develops kinematic discrimination techniques for identifying
**HH → bbτ_h τ_l** (di-Higgs) signal events against tt̄ backgrounds at the
ATLAS detector, using 36.2 fb⁻¹ of ATLAS Open Data (2015–2016 Run 2) and
Monte Carlo simulations.

Di-Higgs production is one of the rarest processes accessible at the LHC.
Observing it would provide a direct measurement of the **Higgs trilinear
self-coupling constant λ₃**, a crucial test of the Standard Model's electroweak
symmetry breaking mechanism. The extremely small production cross-section means
that powerful signal-background discrimination is essential before any
meaningful measurement can be made.

This repository contains ROOT/C++ macros implementing two new penalty functions
— **ΔT2** and **ΔT3** — that extend the tt̄ reconstruction discriminant ψ_T
developed by Cookman & Harris (2019), to handle more complex tt̄ decay
topologies involving tau leptons and additional neutrinos.

---

## Physics Background

### Signal Process

The signal of interest is Standard Model di-Higgs production in the
**bbτ_h τ_l** decay channel:

```
HH → (bb̄)(τ_h τ_l)
         ↓      ↓
      2 b-jets  hadronic tau + light lepton (e/μ) + neutrinos
```

One Higgs decays to a bb̄ pair (detected as two b-tagged jets). The other
decays to two tau leptons — one decaying hadronically (τ_h, detected as a
jet), and one leptonically (τ_l, producing an electron or muon plus neutrinos).
The neutrinos are inferred from the event's missing transverse energy E_T^miss.

### Main Backgrounds

The dominant backgrounds are **dileptonic tt̄ decays**, which produce similar
final states to the signal. This project focuses on three tt̄ decay channels:

| Channel | Decay | Neutrinos |
|---------|-------|-----------|
| **tt1** | tt̄ → bbWW → bbllν_l ν̄_l | 2 |
| **tt2** | tt̄ → bbWW → bblτν_l ν_τ → bbll₂ν_l ν_τ ν_τ₂ ν_l₂ | 4 |
| **tt3** | tt̄ → bbWW → bbτ τ ν_τ^(1) ν_τ^(2) → bbll₂ν_l ν_l₂ ν_τ^(1) ν_τ^(2) ν_τ₂^(1) ν_τ₂^(2) | 6 |

The tt̄ → bbτ_h l channel (tt1 topology) is the single most dominant
background. The presence of up to 6 neutrinos in the more complex topologies
makes full kinematic reconstruction challenging and under-constrained.

### Discrimination Strategy

The overall discriminator is the **ψ penalty function**, first developed by
Langford & Spencer (2017) and extended by Cookman & Harris (2019):

```
Ψ = sqrt(ψ²_bb + ψ²_ττ + ψ²_MET + ψ²_ε + ψ²_ΔR + ψ²_T)
```

Each component penalises an event for inconsistency with the di-Higgs
hypothesis:

- **ψ_bb, ψ_ττ** — invariant masses of bb and ττ pairs vs. Higgs mass (125 GeV)
- **ψ_MET** — consistency of missing transverse energy with neutrino collinearity
- **ψ_ε** — amount of b-jet energy rescaling required
- **ψ_ΔR** — angular separation of b-jets and taus (signal events are boosted
  and collimated; tt̄ events are back-to-back)
- **ψ_T** — how well an event reconstructs under a dileptonic tt̄ hypothesis
  (events consistent with tt̄ are *less* likely to be signal)

Signal events score near zero in ψ; background events score much higher.

---

## This Project's Contributions

Building on the ψ_T discriminant of Cookman & Harris, this project introduces
two new penalty functions to handle tt̄ topologies with more neutrinos:

### ΔT2 — For tt2 events (4 neutrinos)

The two collinear neutrinos from the secondary tau decay (ν_τ₂, ν_l₂) are
approximated as a single composite neutrino whose three-momentum is
proportional to the secondary lepton l₂:

```
p⃗_νt₂νl₂ = β · p⃗_l₂
```

The free parameters minimised are: p⃗_ν_τ, p⃗_ν_l, and the scaling constant β.

### ΔT3 — For tt3 events (6 neutrinos)

Both W bosons decay via tau leptons. Two separate β scaling constants (β₁, β₂)
are used, one for each tau decay chain:

```
p⃗_νl ν_τ₂^(1) = β₁ · p⃗_l
p⃗_νl₂ ν_τ₂^(2) = β₂ · p⃗_l₂
```

The free parameters minimised are: p⃗_ν_τ^(1), p⃗_ν_τ^(2), β₁, and β₂.

### Results

Minimization was performed using **Minuit2** (100 random restarts per event)
on a MC dataset of dileptonic tt̄ events:

| Penalty Function | Events with ΔTᵢ < 400 | Channels Selected |
|-----------------|----------------------|-------------------|
| ΔT (Cookman & Harris) | 71% | tt1 only |
| ΔT2 | 92% | tt1 + tt2 |
| ΔT3 | 94% | tt1 + tt2 + tt3 |

The improvement from ΔT2 and ΔT3 is consistent with the branching ratios of
the respective tt̄ channels (tt1: ~74%, tt2: ~24%, tt3: ~2.5%).

A key empirical finding from minimization experiments: the number of random
restarts dominates convergence quality far more than the choice of mass
distribution shape. Going from 1 → 5 restarts improves convergence by ~180%;
5 → 20 restarts improves it by a further ~30%. Using a Gaussian mass
distribution (W width: ~2 GeV, top width: ~1.3 GeV) rather than fixed masses
improves convergence by only ~7%.

---

## Repository Structure

```
ATLAS-Particle-Physics-Project/
├── delta_t_minimization.C       # ΔT penalty function (tt1 topology, 2 neutrinos)
│                                # Achieves 71% tt̄ event reconstruction efficiency
├── delta_t3_8combinations.C     # ΔT3 penalty function (tt3 topology, 6 neutrinos)
│                                # Iterates over 8 lepton-bjet assignment combinations
├── mini.h                       # Auto-generated ROOT TTree header (not tracked)
└── README.md
```

> **Note:** `mini.h` is auto-generated by ROOT's `MakeClass` from the input
> TTree and is not included in this repository. You will need to generate it
> from your own ROOT file before compiling.

---

## Dependencies

- [ROOT](https://root.cern/) (tested with ROOT 6.x)
- ROOT Math libraries: `Minuit2`, `Math/Minimizer`, `Math/Factory`, `Math/Functor`
- C++11 or later
- Linux environment (developed on Ubuntu with ATLAS Open Data)

---

## How to Run

1. Generate `mini.h` and `mini.C` from your ROOT TTree file:

```bash
root -l your_data_file.root
root [0] T->MakeClass("mini")
```

2. Load and run the macro in ROOT:

```bash
root -l
root [0] .L delta_t3_8combinations.C
root [1] mini t("your_data_file.root", "mini")  // replace with your file
root [2] t.Loop()
```

---

## Background and Prior Work

This project is part of a series of MPhys student projects at the University
of Manchester supervised by Prof. Terry Wyatt:

1. **Langford & Spencer (2017)** — Original b-jet rescaling method and the ψ
   penalty function for HH → bbτ_h τ_l discrimination.
2. **Cookman & Harris (2019)** — Applied ψ to di-Higgs MC events; introduced
   the ΔR discriminant (ψ_ΔR) and the tt̄ reconstruction discriminant (ψ_T)
   via gradient descent; achieved ~70% signal efficiency at 99% background
   rejection using the combined Ψ discriminator.
3. **This project (2019–2020)** — Extended ψ_T with new penalty functions ΔT2
   and ΔT3 to reconstruct the more complex tt̄ → bbτ_h τ_l topology involving
   up to 6 neutrinos, using collinear approximations for secondary tau decay
   neutrino pairs.

---

## Notes

- The code uses fixed W mass (80,379 MeV) and top mass (173,000 MeV).
  Experiments with Gaussian mass distributions were also conducted separately
  (see project report).
- Some files from the project were not version-controlled at the time and are
  no longer available.
- ATLAS Open Data MC samples used include: `WenuWithB`, `WmunuWithB`,
  `VBFH125`, `ZZqqll`, `ZqqZll`.
