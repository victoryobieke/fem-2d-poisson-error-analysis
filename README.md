# fem-2d-poisson-error-analysis
# 2D FEM for Poisson Equation with L² and H¹ Error Analysis

This MATLAB project implements the linear Finite Element Method (FEM) to solve the 2D Poisson problem:

\[
- \Delta u(x, y) = f(x, y), \quad \text{on } \Omega = (0,1)^2, \quad u = 0 \text{ on } \partial\Omega
\]

where the exact solution is known and used for computing the **L² norm** and **H¹ semi-norm** errors over a sequence of mesh refinements.

---

## 🧮 Problem Setup

- **Exact solution:** \( u(x,y) = \sin(\pi x)\sin(\pi y) \)
- **Source term:** \( f(x,y) = 2\pi^2\sin(\pi x)\sin(\pi y) \)
- **Domain:** Unit square \( \Omega = (0,1)^2 \)
- **Boundary conditions:** Homogeneous Dirichlet: \( u = 0 \) on \( \partial\Omega \)

---

## 🔧 Features

- Mesh refinement for `n = [8, 16, 32, 64, 128]`
- Triangular mesh generation via `write_coordinates_dat(n)`
- FEM solver using assembled matrices from `fem2d.m`
- Error computation using 3-point Gaussian quadrature on triangles
- Computes:
  - **L² error**: ‖u - uₕ‖₂
  - **H¹ error**: ‖∇u - ∇uₕ‖
- Automatic convergence rate calculation
- Log-log plot of error vs mesh size \( h = 1/n \)

---

## 📁 Files

- `main_hw3.m`: Main driver script that:
  - Generates meshes
  - Solves the 2D Poisson problem
  - Computes errors and plots convergence
- `error_L2_H1.m`: Function to compute L² and H¹ errors using Gaussian quadrature
- `fem2d.m`: Assembles and solves the FEM linear system
- `write_coordinates_dat.m`: Writes mesh node and connectivity files for each `n`
- `coordinates.dat`, `elements3.dat`: Mesh files used by the solver (generated automatically)

---

## 📈 Output

- **Numerical error table** with L² and H¹ errors and convergence rates
- **Log-log plots** showing \( \mathcal{O}(h^2) \) and \( \mathcal{O}(h) \) reference lines
- Visual diagnostics of FEM accuracy vs mesh resolution

---

## 📦 Requirements

- MATLAB R2018b or later
- No additional toolboxes required

---

## 🙋‍♂️ Author

Victory Obieke  
Oregon State University  
Homework 3 — Finite Element Methods (MTH 656)
