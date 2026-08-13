# Crossflow Model Analysis

MATLAB scripts for analyzing a kinematic crossflow model with base gyre flow + oscillating jet perturbation.

## Prerequisites

- MATLAB
- Custom colormap: [Download here](https://www.mathworks.com/matlabcentral/fileexchange/69470)

## Model Description

Kinematic 2D system based on double-gyre flow (Shadden et al., 2005) with oscillating jet:

$$\dot{\mathbf{x}} = \tilde{\mathbf{u}}(\mathbf{x}, t) + \mathbf{u}'(\mathbf{x}, t)$$

**Base Gyre Flow** $\tilde{\mathbf{u}}$:
```matlab
u = U - π·A·sin(π·f)·cos(π·y)
v = π·A·cos(π·f)·sin(π·y)·df
```
where $f = a·x^2 + b·x$, $a = γ·\sin(ω_g·t)$, $b = 1-2a$

**Oscillating Jet** $\mathbf{u}'$:
```matlab
v' = sech²(b_t(x-a_t))·ε·sin(ω_c·t)
```

**Parameters:**
- Base flow: A=0.1, U=0.2, γ=0.25, ω_g=0 (steady)
- Jet: ε=0.03, ω_c=2π/5 (T_p=5s), b_t=4, a_t=1.4
- Numerics: dt=0.025, dx=0.006, T=15s, domain=[0,2]×[0,1]
- Phase 100: t₀=2.5s (half period)

## Key Concepts

**Phase 100 = T_p/2:**
- Half oscillation period (t₀ = 2.5s)
- Maximum jet effect
- Critical for transport barriers

**FTLE:**
- $\mathrm{FTLE}_{t_0}^{t_0+T} = \frac{1}{T} \log\sqrt{\lambda_{\max}(\nabla\boldsymbol{\Phi}^T \nabla\boldsymbol{\Phi})}$

**MTU Framework:**
- Upper bound: $|\mathbf{x}^\epsilon - \mathbf{x}^0|^2 \leq \Delta_\infty^2 \cdot \mathrm{MS}$
- $\Delta_\infty$: Max jet amplitude $\max_{s \in [t_0, t]} |\mathbf{u}'(\mathbf{x}^0(s), s)|$
- $\mathrm{MS}$: Model sensitivity $\left(\int_{t_0}^t \sqrt{\lambda_s^t} \, ds \right)^2$
- $\zeta$: FTLE perturbation $\frac{1}{T}\left(\log(\Delta_\infty) + \log\left(\frac{\mathrm{MS}}{\lambda_{t_0}^t}\right)\right)$

## Directory Structure

```
crossflow_model/
├── models/
│   ├── crossflowVEC.m               # Combined velocity field
│   ├── gVEC.m                       # Jet perturbation only
│   └── rk4singlestep.m              # RK4 integrator
├── data/                             # Output files
├── simulateCrossflowFTLE.m          # Compute FTLE
├── simulateCrossflowBoundEstField.m # Compute MTU field
├── simulateCrossflowBoundEstPoint.m # Point-wise MTU + validation
└── README.md
```

## Workflow

```matlab
cd crossflow_model

% Step 1: FTLE at phase 100 → FTLE_crossflow_phase100.mat (~1-2 min)
simulate CrossflowFTLE

% Step 2: MTU field → MTU_crossflow_phase100.mat (~5-10 min)
simulateCrossflowBoundEstField

% Step 3: Point-wise MTU + validation (~15-30 min)
simulateCrossflowBoundEstPoint
```

**Total time:** ~20-40 minutes

## What Each Script Does

| Script | Description | Output |
|--------|-------------|--------|
| `simulateCrossflowFTLE.m` | Forward FTLE at phase 100 (T_p/2) | `FTLE_crossflow_phase100.mat` |
| `simulateCrossflowBoundEstField.m` | MTU field (base flow vs jet) | `MTU_crossflow_phase100.mat` |
| `simulateCrossflowBoundEstPoint.m` | Point-wise MTU at 3 locations + validation | `crossflow_point_bounds.mat` |

**Probe Locations:**
- Location 1: (0.48, 0.63)
- Location 2: (1.50, 0.75)
- Location 3: (1.13, 0.918)

## Expected Results

**FTLE Field:**
- Horseshoe structures at x ≈ 0.5, 1.5
- Sharp ridges (material barriers)
- Lobes from jet oscillation
- Typical values: 0-0.3

**MTU Field:**
- High values where jet affects LCS most
- Peaks correlate with FTLE ridges
- Shows jet influence on manifolds

**Point-Wise Validation:**
- Solid lines: Upper bound estimates
- Dotted lines: Actual trajectory errors
- Demonstrates MTU accurately predicts uncertainty

## Citation

**Paper:**
```
Jones, M. R., Klewicki, C. J., Khan, O., Brunton, S. L., & Luhar, M.
"Capturing multiscale interactions in fluid flow via Lagrangian coherent
structures and modal analysis"
```

**Base Model:**
- Shadden et al. (2005). *Physica D*, 212(3-4), 271-304

## Troubleshooting

| Issue | Solution |
|-------|----------|
| `Undefined function 'customcolormap'` | Download from [MATLAB File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/69470) |
| `Cannot find models/` | Check `models/` contains: `crossflowVEC.m`, `gVEC.m`, `rk4singlestep.m` |

