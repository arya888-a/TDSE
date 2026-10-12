# One-Dimensional Time-Dependent Schrödinger Equation (TDSE) Simulation

This repository provides a numerical implementation in Python for solving the Time-Dependent Schrödinger Equation (TDSE) in one spatial dimension. The simulation models the behavior of a Gaussian wave packet as it interacts with a square potential barrier.

## Physical Background

The code numerically solves the dimensionless TDSE:

$$
i \partial_t \Psi(x, t) = \hat{H} \Psi(x, t)
$$

where the Hamiltonian operator $\hat{H}$ is defined as:

$$
\hat{H} = -\frac{1}{2} \partial_x^2 + V(x)
$$

* $\Psi(x, t)$: The wave function of the system.
* $V(x)$: The external potential (in this case, a square barrier).
* $\partial_x^2$: The spatial second derivative computed using finite difference methods.

## Key Features & Components

1. **Spatial and Temporal Grid ($N, L, dx, dt$):** Sets up the spatial grid points within $[-L/2, L/2]$ and defines spatial and temporal steps.
2. **Potential and Initial State:** Constructs a square potential barrier with height $V_0$ and width, and initializes a Gaussian wave packet characterized by initial momentum $p_0$, initial position $x_0$, and width $\sigma$.
3. **Hamiltonian Operator (`apply_H`):** Computes the spatial second derivative using `np.roll` and enforces hard-wall boundary conditions.
4. **Time-Evolution Loop:** Updates the wave function iteratively over time using the forward Euler method.

## Prerequisites

To run this code, you will need the following Python libraries:

* `numpy`
* `matplotlib`
* `IPython` (required for rendering animations in Jupyter environments)

## How to Run

1. Download the main.ipynb and open the code in a suitable environment such as Jupyter Notebook or JupyterLab.
2. Adjust the simulation parameters (such as the number of grid points $N$, box length $L$, and time step $dt$) as needed.
3. Execute the cells sequentially to observe the output animation and the evolution of the wave function norm.

## Technical Notes & Limitations

* **Stability:** The forward Euler method is conditionally stable for quantum differential equations. To prevent divergence in the norm of the wave function, the time step $dt$ must be chosen sufficiently small.
* **Boundary Conditions:** Hard-wall boundary conditions are implemented to contain the wave function within the grid boundaries.
