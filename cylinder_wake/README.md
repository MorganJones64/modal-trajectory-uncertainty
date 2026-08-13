# Cylinder Wake Flow Analysis

MATLAB scripts for analyzing cylinder wake flow at Re=100 using POD, FTLE, and modal-trajectory uncertainty (MTU).

## Prerequisites

- MATLAB (with interpolation toolbox)
- Numerical data from Zenodo
- Custom colormap: [Download here](https://www.mathworks.com/matlabcentral/fileexchange/69470)

## Key Concepts

**Modal Representation:**
- Baseline: $\tilde{u} = \bar{u} + u_1 + u_2$ (mean + first two mode pairs)
- Perturbation: $u' = u_3$ (third mode pair)

**MTU Framework:**
- Quantifies how perturbation modes affect Lagrangian coherent structures
- Upper bound: $|\mathbf{x}^\epsilon - \mathbf{x}^0|^2 \leq \Delta_\infty^2 \cdot \mathrm{MS}$
- $\Delta_\infty$: Maximum perturbation amplitude
- $\mathrm{MS}$: Model sensitivity (integrated stretching)

## Setup

### Download Flow Data

1. Visit: https://zenodo.org/records/5039610
2. Download: `fixed_cylinder_atRe100`
3. Place in: `cylinder_wake/data/`

### Directory Structure

```
cylinder_wake/
├── data/
│   └── fixed_cylinder_atRe100  <- Downloaded file
├── functions/
│   ├── POD.m
│   ├── FTLE.m
│   ├── MTU.m
│   └── MTUbound.m
├── customcolormap/
├── get_dataset.m
├── PODcyl.m
├── calcFTLE.m
├── calcMTU.m
└── calcMTUbound.m
```

## Workflow

```matlab
cd cylinder_wake

% Step 1: Process raw data → cylindervelocity_zen.mat
get_dataset

% Step 2: POD analysis → cylinderPODmodesPair.mat
PODcyl

% Step 3: FTLE fields → bFTLE_*.mat 
calcFTLE

% Step 4: MTU field → MTU_cylinder_wake.mat
calcMTU
NOTE: time to compute may be extensive

% Step 5: Point-wise MTU + Figure 7
calcMTUbound
NOTE: time to compute may be extensive
```

## What Each Script Does

| Script | Description | Output |
|--------|-------------|--------|
| `get_dataset.m` | Reads binary file, interpolates to regular grid | `cylindervelocity_zen.mat` |
| `PODcyl.m` | POD decomposition, extracts mode pairs | `cylinderPODmodesPair.mat` |
| `calcFTLE.m` | Backward FTLE for original & modal flows | `bFTLE_original.mat`, `bFTLE_modalrep.mat` |
| `calcMTU.m` | MTU field ($\tilde{u}=\bar{u}+u_1+u_2$, $u'=u_3$) | `MTU_cylinder_wake.mat` |
| `calcMTUbound.m` | Point-wise MTU + Figure 7 | Bound files + figure |

## Citation

**Paper:**
```
Jones, M. R., Klewicki, C. J., Khan, O., Brunton, S. L., & Luhar, M.
"Capturing multiscale interactions in fluid flow via Lagrangian coherent 
structures and modal analysis"
```

**Citation for Dataset:**
```
Boudina, M. (2021). Numerical simulation data of a two-dimensional flow around 
a fixed circular cylinder [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.5039610
```