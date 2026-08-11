# Crossflow Model Analysis

This directory contains MATLAB scripts for analyzing a kinematic model with a base gyre flow and an oscillating crossflow jet perturbation.

## Model Description

The crossflow model is a 2D kinematic system based on the double-gyre flow from Shadden et al. (2005) with an added oscillating jet perturbation:

$$\dot{\mathbf{x}} = \tilde{\mathbf{u}}(\mathbf{x}, t) + \mathbf{u}'(\mathbf{x}, t)$$
where:
- **$\tilde{\mathbf{u}}(\mathbf{x}, t)$**: Base gyre flow (modal representation)
- **$\mathbf{u}'(\mathbf{x}, t)$**: Oscillating crossflow jet (perturbation)

### Base Gyre Flow Parameters

- **A = 0.1**: Gyre amplitude
- **U = 0.2**: Background flow velocity
- **γ = 0.25**: Asymmetry parameter
- **ω_g = 0**: Gyre oscillation frequency (steady gyres)

### Crossflow Jet Parameters

- **ε = 0.03**: Jet amplitude
- **ω_c = 2π/5**: Jet oscillation frequency (period T_p = 5 seconds)
- **b_t = 4**: Jet width parameter
- **a_t = 1.4**: Jet center position

The jet provides a cross-stream perturbation that oscillates with period T_p = 5 seconds.

## Directory Structure

```
crossflow_model/
├── models/
│   ├── crossflowVEC.m       # Combined flow field (base + perturbation)
│   ├── gVEC.m               # Crossflow jet perturbation
│   └── rk4singlestep.m      # RK4 integrator
├── data/                     # Output data files
├── simulateCurrentftlePUB.m # Compute FTLE at phase 100
├── simulateCurrentmsPUB.m   # Compute MTU (full-field)
├── simulateCurrentbound.m   # Compute point-wise MTU bounds
└── README.md                # This file
```

## Prerequisites

- MATLAB (with basic toolboxes)
- Custom colormap function (`customcolormap`)

## Running the Code

### Step 1: Compute FTLE at Half Period

Compute the forward FTLE field at phase 100 (half oscillation period):

```matlab
cd crossflow_model
simulateCurrentftlePUB
```

**What this does:**
- Computes FTLE starting at t₀ = 2.5 seconds (phase 100 = T_p/2)
- Integration time: T = 15 seconds
- Forward integration to reveal repelling manifolds
- Saves FTLE field to `FTLE_crossflow_phase100.mat`

**Output:**
- **Figure**: Forward FTLE field showing Lagrangian coherent structures
- **File**: `FTLE_crossflow_phase100.mat` containing FTLE field and parameters

**Expected Processing Time:** ~1-2 minutes

### Step 2: Compute Modal-Trajectory Uncertainty

Compute MTU at phase 100:

```matlab
simulateCurrentmsPUB
```

**What this does:**
- Computes MTU at phase 100 (t₀ = 2.5 seconds)
- Tracks particle trajectories with base flow only (no jet)
- Computes perturbation amplitude at intermediate times
- Calculates model sensitivity and trajectory uncertainty bounds

**Output:**
- **Figure**: Scaled MTU field
- **File**: `data/MTU_crossflow_phase100.mat`

**Expected Processing Time:** ~5-10 minutes

**Key Metrics:**
- **Model Sensitivity (MS)**: Integrated stretching from t₀ to t
- **Δ∞**: Maximum jet amplitude along trajectories
- **Scaled MTU**: (MS × Δ∞)² - upper bound on trajectory uncertainty
- **ζ**: FTLE perturbation showing how jet affects FTLE field

### Step 3: Compute Point-Wise MTU Bounds

Compute MTU bounds at three probe locations and validate predictions:

```matlab
simulateCurrentbound
```

**What this does:**
- Computes MTU bounds at 3 locations:
  - Location 1: (0.48, 0.63)
  - Location 2: (1.50, 0.75)
  - Location 3: (1.13, 0.918)
- For each location:
  - Part 1: Computes MTU upper bounds (base flow)
  - Part 2: Computes trajectories WITH jet perturbation
  - Part 3: Computes trajectories WITHOUT jet perturbation
- Plots all results in one figure

**Output:**
- **Figure**: Upper bounds (solid) vs actual errors (dotted) for all 3 locations
- **File**: `data/crossflow_point_bounds.mat`

**Expected Processing Time:** ~15-30 minutes (3 locations)

**Validation:**
- Demonstrates MTU framework accuracy
- Shows upper bounds track actual trajectory errors
- Color-coded by location (green, blue, red)

## Complete Workflow

Run all scripts in sequence:

```matlab
% Navigate to directory
cd crossflow_model

% Step 1: Compute FTLE at phase 100
simulateCurrentftlePUB

% Step 2: Compute full-field MTU at phase 100
simulateCurrentmsPUB

% Step 3: Compute point-wise MTU bounds (optional, for validation)
simulateCurrentbound
```

**Total Time:** 
- Steps 1-2: ~6-12 minutes
- Step 3 (optional): +15-30 minutes

## Key Features

### Phase 100 Analysis

Phase 100 corresponds to t₀ = 100 × dt = 100 × 0.025 = 2.5 seconds, which is:
- **Half the oscillation period**: T_p/2 = 5/2 = 2.5 seconds
- **Maximum jet velocity**: Jet is at peak strength at this phase
- **Critical for LCS**: This phase reveals important transport barriers

### Integration Parameters

- **Time step (dt)**: 0.025 seconds
- **Spatial resolution (dx)**: 0.006
- **Integration duration (T)**: 15 seconds (3 oscillation periods)
- **Domain**: x ∈ [0, 2], y ∈ [0, 1]
- **Grid size**: ~334 × 167 points

### FTLE Computation

The FTLE field is computed as:
$$\mathrm{FTLE}_{t_0}^{t_0+T} = \frac{1}{T} \log\sqrt{\lambda_{\max}(\mathbf{C})}$$

where:
- $\mathbf{C} = \nabla\boldsymbol{\Phi}^T \nabla\boldsymbol{\Phi}$ is the Cauchy-Green strain tensor
- $\nabla\boldsymbol{\Phi}$ is the deformation gradient tensor
- $\lambda_{\max}$ is the maximum eigenvalue

### MTU Computation

The modal-trajectory uncertainty framework computes:

1. **Model Sensitivity (MS)**:
   $$\mathrm{MS}_{t_0}^t = \left(\int_{t_0}^t \sqrt{\lambda_s^t(\mathbf{x}^0(s))} \, ds \right)^2$$
   
   where $\lambda_s^t$ is the maximum eigenvalue at intermediate time $s$

2. **Maximum Perturbation Amplitude**:
   $$\Delta_\infty = \max_{s \in [t_0, t]} |\mathbf{u}'(\mathbf{x}^0(s), s)|$$
   
   where $\mathbf{u}'$ is the jet perturbation field

3. **Trajectory Uncertainty Upper Bound**:
   $$|\mathbf{x}^\epsilon - \mathbf{x}^0|^2 \leq \mathrm{MS}_{t_0}^t \cdot \Delta_\infty^2$$

4. **FTLE Perturbation**:
   $$\zeta_{t_0}^t = \frac{1}{T}\left(\log(\Delta_\infty) + \log\left(\frac{\mathrm{MS}_{t_0}^t}{\lambda_{t_0}^t}\right)\right)$$
   
   This shows how the jet perturbs the FTLE field

## Model Functions

### `crossflowVEC.m`

Computes the combined velocity field:
```matlab
dy = crossflowVEC(t, y, U, A, gamma, omega_g, eps, omega_c, bt, at)
```

**Inputs:**
- `t`: Time
- `y`: Position vector [x; y]
- `U, A, gamma, omega_g`: Base gyre parameters
- `eps, omega_c, bt, at`: Jet perturbation parameters

**Output:**
- `dy`: Velocity vector [dx/dt; dy/dt]

### `gVEC.m`

Computes the crossflow jet perturbation:
```matlab
g = gVEC(t, y, omega, bt, at)
```

Returns the oscillating jet velocity field.

### `rk4singlestep.m`

Fourth-order Runge-Kutta integrator:
```matlab
yout = rk4singlestep(f, dt, t, yin)
```

## Physical Interpretation

### Base Gyre Flow

The base flow consists of two counter-rotating gyres that create:
- **Stable manifolds**: Attracting material lines (backward FTLE ridges)
- **Unstable manifolds**: Repelling material lines (forward FTLE ridges)
- **Hyperbolic points**: Intersections of stable/unstable manifolds

### Crossflow Jet Perturbation

The oscillating jet:
- Creates time-periodic perturbations to the flow
- Generates additional lobes in the manifold structure
- Enhances mixing between gyre regions
- Demonstrates modal-trajectory uncertainty

### Half-Period Significance

At phase 100 (t₀ = T_p/2):
- Jet velocity is at maximum (or minimum, depending on phase)
- Manifold structure shows strongest perturbation effects
- Critical for understanding transport between gyres
- Useful for comparing perturbed vs unperturbed dynamics

## Expected Results

### FTLE Field Features

- **Horseshoe structures**: Around x ≈ 0.5 and x ≈ 1.5
- **Ridge lines**: Sharp FTLE ridges indicating material barriers
- **Lobes**: Additional structures due to jet oscillation
- **Color range**: FTLE values typically 0-0.3

### MTU Field Features

- **High MTU regions**: Where jet most strongly affects particle trajectories
- **Correlation with FTLE**: MTU peaks near FTLE ridges
- **Lobes**: Show where jet creates additional manifold structure
- **ζ field**: Highlights jet's influence on material barriers

### Comparison with Base Flow

When compared to base gyre flow without perturbation:
- **Perturbed manifolds**: Additional corrugation due to jet
- **Enhanced mixing**: Lobes facilitate cross-gyre transport
- **Time-dependence**: Structure varies with initial phase
- **MTU quantification**: Provides upper bound on trajectory deviations

## Applications

This model is useful for:
1. **Understanding mixing mechanisms**: How periodic perturbations affect transport
2. **Testing MTU framework**: Kinematic model with known perturbation
3. **Visualizing LCS**: Clear demonstration of stable/unstable manifolds
4. **Modal analysis**: Separating base flow from perturbation effects
5. **Validating theory**: Known analytical form allows verification of MTU predictions
6. **Educational tool**: Simple system demonstrating mode-manifold interactions

## References

**Base Model:**
- Shadden, S. C., Lekien, F., & Marsden, J. E. (2005). "Definition and properties of Lagrangian coherent structures from finite-time Lyapunov exponents in two-dimensional aperiodic flows." *Physica D*, 212(3-4), 271-304.

**Modal-Trajectory Uncertainty:**
- Kaszás, B., Feudel, U., & Tél, T. (2020). "Leaking in history space: A way to analyze systems subjected to arbitrary driving." *Chaos*, 30(3), 033109.

**Paper (this work):**
- Jones, M. R., Klewicki, C. J., Khan, O., Brunton, S. L., & Luhar, M. "Capturing multiscale interactions in fluid flow via Lagrangian coherent structures and modal analysis"

## Notes

- **Phase selection**: Phase 100 chosen for maximum jet effect at half period
- **Grid resolution**: dx = 0.006 provides good balance between accuracy and speed
- **Integration time**: T = 15 seconds (3 periods) ensures well-developed structures
- **Colormap**: Custom red (forward) or blue (backward) colormap for visualization

## Troubleshooting

### Error: "Undefined function 'customcolormap'"
- Download from: https://www.mathworks.com/matlabcentral/fileexchange/69470-custom-colormap
- Add to path or place in working directory

### Error: "Cannot find models/"
- Ensure `models/` subdirectory exists with required functions
- Check that `crossflowVEC.m`, `gVEC.m`, and `rk4singlestep.m` are present

### Slow computation
- Reduce grid resolution (increase dx)
- Reduce integration time (decrease T)
- Use simpler integration method (Euler instead of RK4)

### FTLE values seem wrong
- Check that phase = 100 (not 1 or other value)
- Verify integration direction (forward vs backward)
- Ensure parameters match expected values

### MTU computation errors
- Ensure FTLE was computed first (creates reference)
- Check that `data/` directory exists
- Verify sufficient memory for large arrays
- If too slow, reduce dx or T

### Point-wise bound computation very slow
- This is normal - computes full trajectories at 3 locations
- Expected time: 15-30 minutes for 3 locations
- Can comment out locations in script to compute fewer points
- Consider running overnight if multiple locations needed

## Data Files

After running all scripts, you will have:

1. **`FTLE_crossflow_phase100.mat`** (from Step 1):
   - `sigma`: FTLE field
   - `x0`, `y0`: Grid coordinates
   - Flow and integration parameters

2. **`data/MTU_crossflow_phase100.mat`** (from Step 2):
   - `scaledMTU`: Upper bound trajectory uncertainty field
   - `zeta`: FTLE perturbation field
   - `mseIntegral`: Model sensitivity
   - `deltaInfty`: Maximum perturbation amplitude
   - `sigma_ftle`: FTLE for baseline system
   - All model parameters

3. **`data/crossflow_point_bounds.mat`** (from Step 3):
   - `all_mseIntegral`: MS at each location (cell array)
   - `all_deltaInfty`: Δ∞ at each location
   - `all_xT_pert`: Particle positions with jet
   - `all_yT_pert`: Particle positions with jet
   - `all_xT_clean`: Particle positions without jet
   - `all_yT_clean`: Particle positions without jet
   - `locations`: Probe point coordinates
   - `locationNames`: Names of each location

## Contact

For questions about this code, please open an issue in the repository.
