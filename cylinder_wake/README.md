# Cylinder Wake Flow Analysis

This directory contains MATLAB scripts for analyzing flow data around a cylinder at Re=100.

## Prerequisites

- MATLAB (with interpolation toolbox)
- Flow data from Zenodo

## Setup Instructions

### 1. Download the Flow Data

Before running any scripts, you must download the required flow dataset:

1. Visit the Zenodo repository: https://zenodo.org/records/5039610
2. Download the file `fixed_cylinder_atRe100`
3. Place the downloaded file in the `data/` subdirectory:
   ```
   cylinder_wake/
   ├── data/
   │   └── fixed_cylinder_atRe100  <- Place the file here
   ├── subgrid_dataset.m
   └── README.md
   ```

### 2. Verify Directory Structure

Ensure your directory structure looks like this:
```
cylinder_wake/
├── data/
│   └── fixed_cylinder_atRe100
├── functions/
│   ├── POD.m
│   ├── CSEsb.m
│   ├── FTLEsb.m
│   └── MS2bound.m
├── customcolormap/          (if available)
├── subgrid_dataset.m
├── PODcyl.m
├── calcCSEsandbox.m
├── calcFTLEsandbox_vid.m
└── ... (other MATLAB scripts)
```

## Running the Code

The complete analysis workflow consists of four main steps. Run them in sequence:

### Step 1: Generate the Subgrid Dataset

The main data processing script is `subgrid_dataset.m`. It performs two operations:

**Part 1: Read Flow Data**
- Reads the binary file `fixed_cylinder_atRe100`
- Extracts Reynolds number, velocity fields (U, V), pressure, and coordinates
- Saves raw data to `flow_data.mat`

**Part 2: Interpolate onto Regular Subgrid**
- Loads the raw flow data
- Interpolates scattered data onto a regular grid
- Extracts a subdomain of interest: x ∈ [-20, 20], y ∈ [-7, 7]
- Saves processed data to `cylindervelocity_zen.mat`

To run:
```matlab
cd cylinder_wake
subgrid_dataset
```

### Output Files

After running `subgrid_dataset.m`, you will have:

1. **`flow_data.mat`** - Raw flow data
   - `Re`, `Ur`: Reynolds number and reduced velocity
   - `times`: Time vector
   - `X`, `Y`: Scattered node coordinates
   - `U`, `V`: Velocity components at scattered nodes

2. **`cylindervelocity_zen.mat`** - Processed gridded data
   - `u`, `v`: Velocity fields on regular grid (ny × nx × nt)
   - `xMat`, `yMat`: 2D coordinate meshgrids
   - `x`, `y`: 1D coordinate vectors
   - `dx`, `dy`, `dt`: Grid spacing
   - `nx`, `ny`, `nt`: Grid dimensions

### Step 2: Perform POD Analysis

After generating the gridded dataset, run POD analysis:

```matlab
PODcyl
```

This performs Proper Orthogonal Decomposition on the velocity data and produces:

1. **`cylinderPOD.mat`** - POD decomposition
   - `Zeta`, `Sigma`, `Xi`: POD matrices

2. **`cylinderPODmodesSingle.mat`** - Individual modes
   - `u11`, `u12`, `u21`, `u22`, etc.: Mode velocities
   - `umean`, `vmean`: Mean flow

3. **`cylinderPODmodesPair.mat`** - Mode pairs  
   - `u1`, `u2`, `u3`, `u4`, `ut`: Mode pair velocities
   - `umean`, `vmean`: Mean flow

### Step 3: Visualize Modal Representation and Perturbation

Visualize the modal representation and perturbation fields (Figure 5 from the paper):

```matlab
plotdynsys
```

This script generates:
- **Figure 1**: Modal representation $\tilde{u} = \bar{u} + u_1 + u_2$ (streamwise velocity)
- **Figure 2**: Perturbation $u' = u_3$ (streamwise velocity)

**Output:**
- Visual comparison between the baseline modal representation and the perturbation mode
- Shows how the third mode pair acts as a perturbation to the dominant wake dynamics

### Step 4: Compute and Visualize FTLE Fields

Compute backward FTLE (Finite-Time Lyapunov Exponent) fields:

```matlab
calcFTLE
```

This script:
1. Computes backward FTLE for the original flow field
2. Computes backward FTLE for the modal representation ($\tilde{u} = \bar{u} + u_1 + u_2$)
3. Visualizes both FTLE fields for comparison

**Output Files:**
1. **`bFTLE_original.mat`** - FTLE for original flow
   - `sigma`: Maximum stretching values
   - `xPos`, `yPos`: Particle positions
   - Integration parameters

2. **`bFTLE_modalrep.mat`** - FTLE for modal representation
   - `sigma_modal`: Maximum stretching values
   - `xPos_modal`, `yPos_modal`: Particle positions
   - Integration parameters

**Figures:**
- **Figure 1**: Backward FTLE field of original flow
- **Figure 2**: Backward FTLE field of modal representation
- Reveals Lagrangian coherent structures (LCS) and von Kármán vortex street

### Step 5: Compute Modal-Trajectory Uncertainty

Compute backward modal-trajectory uncertainty (MTU) to quantify how mode $u_3$ influences the LCS:

```matlab
calcMTU
```

This script:
1. Constructs modal representation: $\tilde{u} = \bar{u} + u_1 + u_2$
2. Defines perturbation: $u' = u_3$
3. Computes MTU metrics including scaled MTU field
4. Visualizes the scaled MTU field

**Output Files:**
- **`MTU_cylinder_wake.mat`** - MTU computation results
  - `sigma_ftle`: FTLE values for baseline system
  - `cseIntegral`: Integrated sensitivity
  - `deltaInfty`: Maximum perturbation amplitude
  - `xPos`, `yPos`: Particle positions
  - Integration parameters

**Figures:**
- **Figure 1**: Scaled modal-trajectory uncertainty field
- Highlights regions where mode $u_3$ most strongly influences LCS
- Provides insight into multi-scale interactions between modes and manifolds

**Note:** For computing MTU with different configurations (forward/backward, multiple frames), use `calcMTU_general.m`

### Step 6: Compute Point-Wise MTU Analysis and Generate Figure 7

Generate complete trajectory uncertainty analysis (Figure 7 in the paper):

```matlab
calcMTUboundSweep
```

This script performs the complete point-wise analysis in four parts:

**Part 1: Compute MTU bounds at probe points**
- Probe points: Selected along the vortex street at different streamwise positions
- Computes $\Delta_{\infty}$ and integrated sensitivity at each point
- Saves `bounddatasave_backwards_p*.mat` files

**Part 2: Compute particle trajectories (original flow)**
- Tracks particles backward in time using full original velocity field
- Saves `particlePos_backwards_p*_original.mat` files

**Part 3: Compute particle trajectories (modal representation)**
- Tracks particles using modal representation: $\tilde{u} = \bar{u} + u_1 + u_2$
- Saves `particlePos_backwards_p*_fmnm1m2.mat` files

**Part 4: Plot upper bounds vs trajectory errors (Figure 7)**
- Loads all computed data
- Compares upper bound estimates to actual trajectory errors
- Shows how well MTU predicts trajectory uncertainty
- Demonstrates correlation between $\zeta$ values and trajectory errors

**Output Files:**
- **`bounddatasave_backwards_p*.mat`** - MTU bound data for each point
  - `cseIntegral`: Integrated sensitivity over time
  - `deltaInfty`: Maximum perturbation amplitude
  - Point-specific parameters

- **`particlePos_backwards_p*_original.mat`** - Particle positions (original flow)
  - `xPos`, `yPos`: Particle trajectories
  - `sigma_ftle`: FTLE values

- **`particlePos_backwards_p*_fmnm1m2.mat`** - Particle positions (modal representation)

**Figure:**
- Trajectory uncertainty upper bounds (solid lines) vs actual errors (dotted lines)
- Shows predictive capability of MTU framework
- Validates model sensitivity approach for fluid flows

---

## Complete Workflow Summary

Run the following commands in sequence:

```matlab
% Navigate to the cylinder_wake directory
cd cylinder_wake

% Step 1: Generate gridded velocity data from raw file
subgrid_dataset

% Step 2: Perform POD analysis
PODcyl

% Step 3: Visualize modal representation and perturbation
plotdynsys

% Step 4: Compute and visualize FTLE fields
calcFTLE

% Step 5: Compute modal-trajectory uncertainty
calcMTU

% Step 6: Complete point-wise MTU analysis and Figure 7
calcMTUboundSweep
```

**Expected Processing Time:**
- Step 1: 10-20 minutes (interpolation)
- Step 2: 2-5 minutes (POD computation)
- Step 3: < 1 minute (visualization)
- Step 4: 20-40 minutes (FTLE computation for both fields)
- Step 5: 30-60 minutes (MTU computation)
- Step 6: 10-20 minutes (point-wise MTU, particle trajectories, and plotting)

**Total:** ~65-140 minutes depending on system

## Dataset Information

- **Source**: Zenodo record 5039610
- **Flow Configuration**: Flow past a fixed cylinder
- **Reynolds Number**: Re = 100
- **Grid Resolution**: 2000 × 2000 (full), extracted to subdomain
- **Subdomain**: x ∈ [-20, 20], y ∈ [-7, 7]

## Analysis Methods

### Proper Orthogonal Decomposition (POD)

POD identifies energetically dominant coherent structures by decomposing the velocity field into orthogonal spatial modes ranked by energy content. For the cylinder wake:
- Mode pairs $u_1$ and $u_2$ capture the primary vortex shedding frequency
- Higher mode pairs ($u_3$, $u_4$, etc.) represent harmonics and smaller-scale structures
- The modal representation $\tilde{u} = \bar{u} + u_1 + u_2$ captures the dominant wake dynamics
- Mode $u_3$ serves as a perturbation to study modal-trajectory uncertainty

### Finite-Time Lyapunov Exponent (FTLE)

FTLE identifies Lagrangian coherent structures (LCS) by tracking fluid particle trajectories:
- **Backward FTLE** reveals attracting manifolds (stable structures)
- **Forward FTLE** reveals repelling manifolds (unstable structures)
- Ridges in the FTLE field correspond to material surfaces that organize particle motion
- For cylinder wake, FTLE reveals the von Kármán vortex street structure

### Modal-Trajectory Uncertainty (MTU)

This analysis (described in the paper) examines how specific POD modes influence LCS:
- Based on the model sensitivity framework by Kaszás et al. (2020)
- Computes the upper-bound trajectory uncertainty due to modal perturbations
- **Scaled MTU field**: Highlights regions where perturbation modes most strongly affect LCS
- **FTLE perturbation ($\zeta$)**: Shows how modes perturb the FTLE field
- Quantifies multi-scale interactions between Eulerian modes and Lagrangian structures

**Key Metrics:**
- $\Delta_{\infty}$: Maximum perturbation amplitude along trajectories
- $\mathrm{MS}$: Model sensitivity (integrated stretching)
- Scaled MTU: $\log(\Delta_{\infty}^2 \mathrm{MS}) / (2|T|)$

For the cylinder wake:
- Baseline system: $\tilde{u} = \bar{u} + u_1 + u_2$
- Perturbation: $u' = u_3$ (third mode pair)
- Reveals how harmonic modes influence vortex street dynamics

## Important Notes

- All scripts display progress messages during execution
- Large computations (interpolation, FTLE) can take significant time
- MATLAB v7.3 format is used for large data files (> 2GB)
- Ensure sufficient RAM (recommended: 16GB+) for FTLE computations
- The `customcolormap/` directory is required for proper visualization

## References and Citation

If you use this code or data, please cite:

**Original Dataset:**
```
Zenodo Record: https://zenodo.org/records/5039610
```

**Paper:**
```
Jones, M. R., Klewicki, C. J., Khan, O., Brunton, S. L., & Luhar, M.
"Capturing multiscale interactions in fluid flow via Lagrangian coherent 
structures and modal analysis"
```

**Key References:**
- POD: Berkooz, G., Holmes, P., & Lumley, J. L. (1993). *Annual Review of Fluid Mechanics*
- FTLE: Shadden, S. C., et al. (2005). *Physica D*
- Modal Sensitivity: Kaszás, B., et al. (2020). *Chaos*

## Troubleshooting

### Error: "Unable to open file 'data/fixed_cylinder_atRe100'"
- Make sure you downloaded the file from Zenodo
- Verify the file is in the `data/` subdirectory
- Check file permissions

### Error: "Undefined function or variable 'FTLE'"
- Ensure `functions/` directory contains `FTLE.m`
- Check that paths are added correctly with `addpath('functions/')`

### Error: "Undefined function 'customcolormap'"
- Ensure the `customcolormap/` directory exists
- Download from: https://www.mathworks.com/matlabcentral/fileexchange/69470-custom-colormap
- Add to path: `addpath('customcolormap/')`

### Memory Issues
- Interpolation and FTLE computations use significant memory
- Recommended: 16GB+ RAM for full analysis
- If memory errors occur:
  - Reduce `Nx` and `Ny` in `subgrid_dataset.m`
  - Reduce `ROInx` and `ROIny` in `calcFTLE.m`
  - Close other applications during computation

### Long Processing Times
- FTLE and MTU computations are computationally intensive
- Integration method `'RK4'` is accurate but slower than `'Euler'`
- MTU requires computing particle trajectories at multiple integration times
- Consider running overnight for full analysis
- Progress messages will indicate computation status

### MTU-Specific Issues
- MTU computation is more expensive than FTLE (requires trajectory tracking at intermediate times)
- Currently only `extrap = true` mode is fully implemented
- For custom configurations (forward/backward, multiple frames), use `calcMTU_general.m`
- Ensure `functions/MTU.m` is in the path

## Contact

For questions about this code, please open an issue in the repository.
-