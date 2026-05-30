# 2D Finite Volume Euler Solver for Airfoils

This project is a 2D Computational Fluid Dynamics (CFD) solver written in C++ (Visual Studio 2022). It solves the compressible Euler equations using a cell-centered Finite Volume Method (FVM). Based on the periodic wraparound logic in the `xi` direction (`wrapXi`) and the output file names, it is specifically designed for structured O-grid meshes around aerodynamic bodies. 

## 1. Physics & Numerical Methods

* **Governing Equations:** 2D Compressible Euler Equations (Inviscid flow).
* **Spatial Discretization:** Cell-centered Finite Volume Method.
* **Convective Flux Scheme:** Roe's Flux-Difference Splitting (`RoeFlux::computeFlux`).
* **Spatial Reconstruction:** Up to 2nd-Order TVD (Total Variation Diminishing).
* **Slope Limiter:** Minmod limiter to prevent spurious oscillations near shocks.
* **Time Integration:** Explicit 3rd-Order TVD Runge-Kutta (Shu-Osher form).
* **Boundary Conditions:**
    * **Solid Wall:** Slip condition utilizing Riemann invariants and isentropic relations.
    * **Far-Field:** Characteristic-based boundary condition using 1D Riemann invariants.

## 2. Core Data Structures

| Structure | Description |
| :--- | :--- |
| `StateVec` | Represents the conservative state vector: `U = [rho, rho*u, rho*v, rho*E]^T`. Contains helper methods to extract primitives (`u, v, p, H, a`) and overloaded operators for vector math. |
| `Field2D<T>` | A flattened 1D `std::vector` wrapper acting as a contiguous 2D array. Accessed via `(i, j)`. Used for states, cell volumes, and face normals. |
| `FaceNormal` | Stores the outward-facing normal components (`nx`, `ny`) and the `length` of a cell face. |
| `Point2D` | Simple struct for raw X and Y node coordinates. |

## 3. Project Architecture

The repository is modularized into the following namespaces/directories:

* **Core (`Type.h`):** Defines global constants (e.g., gamma = 1.4) and the fundamental data types.
* **IO (`MeshLoader`, `TecplotWriter`):** Handles loading nodal coordinates from structured grid files (Plot3D/Text) and exporting cell-centered primitive variables to Tecplot format (`.plt`).
* **Boundaries (`Boundary`):** Computes physical fluxes directly at the boundaries without generating ghost cells.
* **Numerics (`RoeFlux`):** Contains the exact math for Roe averaging, eigenvalue/eigenvector dissipation, and flux assembly.
* **Solver (`EulerSolver2D`):** The orchestrator. Manages the grid variables, handles spatial reconstruction (`xi` and `eta` directions), extrapolates boundary states, and computes the complete spatial residual matrix.
* **TimeIntegration (`RungeKutta3`):** Applies the explicit TVD RK3 step, managing the fractional memory buffers (`U_n`, `dU`) and dynamically fetching the local CFL-restricted time step.

## 4. Execution Pipeline & Data Flow

1. **Initialization:** `MeshLoader` reads the `.plt` mesh, calculates cell volumes using cross-products of diagonals, and computes face normals. `EulerSolver2D` is initialized with free-stream conditions (`U_inf`).
2. **Time-Stepping Loop:** `RungeKutta3::step` is called.
3. **Local Time Step:** Calculates the minimum time step based on the provided CFL number and maximum local wave speeds.
4. **Residual Computation (Per RK Stage):**
    * Zeroes out the residual field.
    * Loops over internal `xi` and `eta` faces.
    * Reconstructs Left (L) and Right (R) states at the face interfaces using Minmod.
    * Passes L and R states to `RoeFlux` to get the interface flux.
    * Adds/subtracts fluxes to adjacent cell residuals.
    * Computes and applies boundary fluxes (Solid Wall and Far-field).
5. **State Update:** Updates `EulerSolver2D::U` using the RK3 fractional steps.
6. **Export:** Periodically calls `TecplotWriter` to dump the updated primitive variables for visualization.

## 5. Context Notes for AI Assistance

* **Grid Topology:** The `xi`-direction (`i` index) utilizes a `wrapXi` function, implying a wrap-around mesh (like an O-grid). The `eta`-direction (`j` index) is bounded by the solid wall at `j=0` and the far-field at `j = N_y - 1`.
* **Memory Layout:** `Field2D` is strictly zero-indexed for cell centers. `nx` and `ny` in the solver represent the number of *cells*, while the mesh loader reads nodes.
* **Flux Addition:** Residuals represent `dU/dt`. Fluxes are added to the right cell and subtracted from the left cell based on the standard finite volume sign convention.