# Time-Dependent Schrödinger Equation Simulation

A numerical simulation of the **1D time-dependent Schrödinger equation (TDSE)** for a Gaussian wave packet interacting with a square potential barrier.

The wave function is discretized in space using a **finite-difference method** and evolved in time using the **Forward Euler method**. The simulation also tracks the wave-function norm to monitor numerical stability and generates an animation of the probability density ( |\Psi(x,t)|^2 ).

### Features

* Gaussian wave-packet initialization and normalization
* Square potential barrier
* Finite-difference approximation of the second spatial derivative
* Forward Euler time evolution
* Hard-wall boundary conditions
* Norm conservation/stability monitoring
* Animated visualization of wave-packet evolution
