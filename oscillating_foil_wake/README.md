# Oscillating Foil Wake Analysis

MATLAB scripts for analyzing PIV measurements of an oscillating NACA-0012 foil wake at $Re_c=11{,}000$ and $St=0.29$ using optimized DMD, finite-time Lyapunov exponents (FTLE), and modal-trajectory uncertainty (MTU).

## Prerequisites

- MATLAB with support for `griddedInterpolant`, `pagemtimes`, and argument validation
- Oscillating hydrofoil wake data from Hugging Face
- The opt-DMD and custom-colormap functions included in this directory

## Modal System

Following equation `osc-sys1` in the paper:

- Baseline: $\tilde{\mathbf{u}}=\bar{\mathbf{u}}+\mathbf{u}_1$
- Perturbation: $\mathbf{u}'=\mathbf{u}_2$
- Perturbed system: $\tilde{\mathbf{u}}+\mathbf{u}'=\bar{\mathbf{u}}+\mathbf{u}_1+\mathbf{u}_2$

The mean and four unsteady mode pairs are obtained with constrained optimized DMD. The eigenvalues are restricted to the imaginary axis so the fitted modes are periodic and have no growth or decay.

## Download Flow Data

1. Visit [morganrj071/oscillating-hydrofoil-wake-Re11000](https://huggingface.co/datasets/morganrj071/oscillating-hydrofoil-wake-Re11000).
2. Download `oscfoil_data.mat` and `foilcoords.mat`.
3. Create `oscillating_foil_wake/data/`.
4. Place both files in that folder.

The expected layout is:

```text oscillating_foil_wake/README.md
oscillating_foil_wake/
├── data/
│   ├── oscfoil_data.mat
│   └── foilcoords.mat
├── functions/
│   ├── FTLEsb.m
│   ├── CSEsb.m
│   └── optdmdsrc/
├── customcolormap/
├── OptDMD.m
├── calcFTLE.m
├── calcMTU.m
└── README.md
```

## Workflow

Run the scripts in order from the `oscillating_foil_wake` directory:

```matlab oscillating_foil_wake/README.md
cd oscillating_foil_wake

% Step 1: Fit opt-DMD and reconstruct the modal velocity fields.
OptDMD

% Step 2: Compute original, baseline, and perturbed FTLE fields.
calcFTLE

% Step 3: Compute the MTU field for baseline umean+u1 and perturbation u2.
% This is the most computationally expensive step.
calcMTU
```

The scripts are intentionally separate because each stage is computationally expensive and produces reusable `.mat` files. A completed stage does not need to be repeated when working on a later stage or changing only plot settings.

## Scripts and Outputs

| Script | Description | Main output |
|---|---|---|
| `OptDMD.m` | Fits nine constrained opt-DMD modes and reconstructs the mean plus four mode pairs | `data/optdmd_coeffs.mat`, `data/opdmd_modes.mat` |
| `calcFTLE.m` | Computes backward FTLE fields for the original, baseline, and perturbed systems | `data/bFTLE_original.mat`, `data/bFTLE_baseline.mat`, `data/bFTLE_perturbed.mat` |
| `calcMTU.m` | Computes the backward MTU field for $\tilde{u}=\bar{u}+u_1$ and $u'=u_2$ | `data/MTU_oscillating_foil.mat` |

## Computation Notes

- One foil period contains approximately 99 samples (`pnt = 99`).
- `calcFTLE.m` uses a backward integration interval of approximately half a period.
- `calcMTU.m` uses a backward integration interval of approximately 1.5 periods.
- Both calculations use fourth-order Runge-Kutta integration.
- The MTU calculation repeats flow-map integrations over intermediate times and can take substantially longer than the FTLE calculation.
- To test the workflow quickly, reduce `ROInx`, `ROIny`, and/or `abs(tLength)` before running the publication-resolution calculation.

## Boundary Approximation

The measured PIV domain is spatially limited. Outside the domain, the scripts use the assumptions described in the paper:

- Baseline streamwise velocity approaches the normalized freestream value $u_\infty=1$.
- Baseline transverse velocity is zero.
- The perturbation velocity is zero.

## MTU Quantities

The saved fields include:

- `sigma_ftle`: square root of the maximum Cauchy-Green eigenvalue for the baseline flow
- `cseIntegral`: integrated modal sensitivity stretching
- `deltaInfty`: maximum perturbation amplitude along baseline trajectories
- `xPos`, `yPos`: advected particle positions

The plotted scaled MTU is

$$
\frac{1}{2|T|}\log\left[\Delta_\infty^2\left(\int\sqrt{\lambda_s^t}\,ds\right)^2\right].
$$

## Citation

```text oscillating_foil_wake/README.md
Jones, M. R., Klewicki, C. J., Khan, O., Brunton, S. L., & Luhar, M.
"Capturing multiscale interactions in fluid flow via Lagrangian coherent
structures and modal analysis"
```

The experimental oscillating-foil dataset is associated with the wake measurements discussed by Jones, Kanso, and Luhar and is distributed at the Hugging Face link above.
