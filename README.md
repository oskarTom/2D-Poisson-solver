# 2D Poisson Solver

## Overview

This repository provides a parallelized C++ implementation for solving the two-dimensional Poisson's equation on a unit square, utilizing the Successive Over-Relaxation (SOR) algorithm. The code is designed to run efficiently on multiple CPU cores using MPI for distributed memory parallelism. 

## Mathematical Background

The Poisson's equation on a unit square is given by:

$$
\frac{\partial^2}{\partial x^2}f(x,y) + \frac{\partial^2}{\partial y^2}f(x,y) = g(x,y), \quad (x, y) \in [0,1]^2
$$

with specified boundary conditions and a known function $g(x, y)$. The solution is discretized, and central finite differences are used to approximate the derivatives.

### Iterative Methods

- **Jacobi Method:** Updates each grid point using the average of its neighbors from the previous iteration.
- **Gauss-Seidel Method:** Utilizes the most recently updated values for faster convergence.
- **Successive Over-Relaxation (SOR):** Introduces a relaxation parameter $\gamma$ to accelerate convergence:

$$
  f_{i,j}^{n+1} = (1-\gamma)f_{i,j}^{n} + \frac{\gamma}{4}(f_{i+1,j}^{n}+ f_{i-1,j}^{n+1}+ f_{i,j+1}^{n}+ f_{i,j-1}^{n+1} ) -\frac{\gamma}{4N^2}g_{i,j}
$$

## Parallelization Approach

The code employs the red-black SOR algorithm to enable parallel updates:

- The grid is split into "red" and "black" cells (checkerboard pattern).
- Each parallel process (1, 2, or 4 supported) manages a subgrid.
- At each iteration:
  1. Update red sublattice.
  2. Synchronize border values with neighboring processes.
  3. Update black sublattice.
  4. Synchronize again.
- Subgrids are merged at the end of computation.

The parallelization uses MPI for process communication and synchronization.

## Usage Instructions

### Requirements

- MPI library (e.g., `OpenMPI` or `MPICH`)
- C++ compiler (`g++` recommended)

### Compilation

```bash
mpic++ -o poisson_solver main.cpp
```


### Running the Solver

The program can be run using 1, 2, or 4 processes:

```bash
mpirun -np 4 ./poisson_solver
```

**Note:** The code is currently configured to work only with 1, 2, or 4 processes.

### Customizing the Problem

- The function $g_{i,j}$ (right-hand side of Poisson's equation) and boundary conditions are hardcoded. To solve for different problems, modify the relevant functions in the code.
- The grid size $N$ and the number of iterations can be set in the source code.


## Code Structure

Key functions:

- `double f(double f, double l, double r, double d, double u, double g, int N)`: Updates a grid point using SOR.
- `void red_black(double** M, int N, int iters, int id, int ntasks)`: Main parallelized SOR routine.
- MPI is initialized and the grid is set up in `main`.

Processes communicate border values for correct parallel updates and finally merge results.

## Benchmarking

The code can be used to simulate physical scenarios such as the electrostatic potential due to a point charge, providing a comparison with the analytical Green's function. While the solver converges well near the source, the value at the source may not match the analytical value exactly due to discretization.

## Limitations and Possible Improvements

- The number of processes is limited to 1, 2, or 4.
- The right-hand side $g_{i,j}$ and boundary conditions are hardcoded.
- No automated benchmarking or visualization is included.
- Merging and communication routines could be optimized.
- The code could be refactored for generality and modularity.

## References

- William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery, *Numerical Recipes in C*, Second Edition (1992), Cambridge University Press.

---

For further mathematical and technical details, see `report.tex`.
