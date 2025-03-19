# Closed Systems Module

## Overview

In this section, we outline the main functionality of the `closed_systems` module. 

First, we will provide some conceptual explanations that provide the user with a necessary background to understand the code. Then we provide some illustrative examples that demonstrate how the code can be used. Finally, we provide the source code as an API reference to the source code.

## Examples and Introductory Concepts 

Before we look at doing Quantum Dynamics on a quantum computer, we'll start out by looking at some ways that we can do quantum dynamics on a classical computer! This is an important step as it will familiarize you with the general ingredients of a quantum dynamics simulation and will also provide us with a means of validating the results obtained from a quantum computer.

### **Your first simulation: The Quantum Harmonic Oscillator**

#### Propagation in the Ladder Basis

We begin by showing a simple example of how to compute dynamics of an initial state using QuTiP's mesolve() function. When we run dynamics, we must do the following:

1. **Define the initial state** $\left| \alpha \right\rangle$. In this example, our initial state is defined as a coherent state with amplitude $\alpha = x_{0} + i p_{0}$, which can be expressed in the Fock Basis as:
     
    $$ \left| \alpha \right\rangle = \frac{\alpha^{n}}{\sqrt{n!}} e^{-\frac{1}{2}\left| \alpha \right|^{2}} \left| n \right\rangle $$
     
2. **Define the Hamiltonian** $H$. In this example, our Hamiltonian is the familiar quantum harmonic oscillator Hamiltonian, defined in terms of creation and annihilation operators as:

    $$ H =  \hbar \omega \left( \hat{a}^{\dagger} \hat{a} + \frac{1}{2} \right) $$

3. **Define the propagation time step $t$ and the number of time steps $n$ for which to compute the wavefunction.**

4. **Compute the time-evolved wavefunction at each step as:**

    $$ \left| \alpha (t_{i+1}) \right\rangle = e^{-\frac{i}{\hbar} H t} \left| \alpha(t_{i}) \right\rangle $$

All this can be done using QFlux as follows: 

```python
# Import the package and relevant modules
import qflux
from qflux.closed_systems import DynamicsCS
# Instantiate our Closed-Systems Dynamics Class
qho_dyn_obj = DynamicsCS(n_basis=128, xo=1.0, po=0.0, mass=1.0, omega=1.0)
# Define our coordinate x and p operators
qho_dyn_obj.set_coordinate_operators()
# Initialize the ladder operators
qho_dyn_obj.initialize_operators()
# Define the default initial state (note that custom initialization is also supported)
qho_dyn_obj.set_initial_state()
# Define some parameters for the time evolution
total_time = 20.0
N_steps = 400
qho_dyn_obj.set_propagation_time(total_time, N_steps)
# Set the Potential/Hamiltonian for our object, in this case using the pre-defined 'harmonic' oscillator potential
qho_dyn_obj.set_hamiltonian(potential_type='harmonic')
# Propagate using the QuTiP sesolve method
qho_dyn_obj.propagate_qt()
```

We can validate our result by computing and plotting the expectation values of $x$ and $p$ as a function of time and comparing to the analytic results: 

```python
import qutip as qt
import numpy as np
import matplotlib.pyplot as plt 

# Compute expectation values <x> and <p>
exp_x_qt = qt.expect(qho_dyn_obj.x_op, qho_dyn_obj.dynamics_results_op.states)
exp_p_qt = qt.expect(qho_dyn_obj.p_op, qho_dyn_obj.dynamics_results_op.states)

exp_x_ana = [ qho_dyn_obj.xo*np.cos(qho_dyn_obj.omega*t) + (qho_dyn_obj.po/qho_dyn_obj.mass/qho_dyn_obj.omega)*np.sin(qho_dyn_obj.omega*t) for t in qho_dyn_obj.tlist]
exp_p_ana = [ qho_dyn_obj.po*np.cos(qho_dyn_obj.omega*t) -qho_dyn_obj.xo*qho_dyn_obj.omega*qho_dyn_obj.mass*np.sin(qho_dyn_obj.omega*t)  for t in qho_dyn_obj.tlist]

# Plot the final result
plt.figure(figsize=(9, 6.5))
plt.plot(qho_dyn_obj.tlist, exp_x_qt, label=r'$\left\langle x \right\rangle$ (QuTiP)', color='dodgerblue')
plt.plot(qho_dyn_obj.tlist, exp_x_ana, label=r'$\left\langle x \right\rangle$ (Analytic)',
         lw=0, marker='x', markevery=10, color='dodgerblue', ms=8)
plt.plot(qho_dyn_obj.tlist, exp_p_qt, label=r'$\left\langle p \right\rangle$ (QuTiP)', color='crimson')
plt.plot(qho_dyn_obj.tlist, exp_p_ana, label=r'$\left\langle p \right\rangle$ (Analytic)',
         lw=0, marker='x', markevery=10, color='crimson', ms=8)
# plt.plot(qho_dyn_obj.tlist, exp_p_grid, color='orange')

plt.ylim(-1.55, 1.55)

plt.legend(ncols=2, loc='upper center')
plt.hlines([-1, 0, 1], min(qho_dyn_obj.tlist), max(qho_dyn_obj.tlist), ls='--', lw=0.85, color='tab:grey', zorder=2)
```

![expectation_vals_ho](images/QHO-Expectation-Values.png)

Hopefully the agreement in this plot convinces you that we're doing something correct! 

#### Propagation in the Coordinate Basis

We can also propagate a wavefunction in a coordinate-grid representation using the so-called Split-Operator Fourier Transform (SOFT) method. 

$\renewcommand{\intertext}[1]{\\\ \textrm{#1}\\}$

Here, we compute the time-evolution of a wavepacket defined in the position basis ($\psi(x)$) according to the Split-Operator Fourier Transform (SOFT) method. This differs slightly from the approach utilized in the last section in the following ways:

- As mentioned above, we will describe the wavefunction in terms of the position. To do this, we must define a closed range of positions $x$ and momenta $p$ and discretize over some finite number of points (analagous to the finite number of Fock states considered above).
- We will compute the time evolution as:

$$ \left| \psi (t) \right\rangle =  e^{- \frac{i}{\hbar} {H} t} \left| \psi(0) \right\rangle $$

Writing $H$ in terms of Kinetic and Potential energy $H = T + V$:

\begin{align*}
        \left| \psi (t) \right\rangle &\approx \lim\limits_{N\to\infty} \left[e^{\frac{-ip^2t}{2m\hbar N}}\ e^{\frac{-i{V}({x})t}{\hbar N}}\right]^N \space |\psi(0)\rangle \\
                                  &\approx\lim\limits_{N\to\infty} \left[e^{\frac{-i{V}({x})t}{2\hbar N}}\space e^{\frac{-ip^2t}{2m\hbar N}}\space e^{\frac{-i{V}({x})t}{2\hbar N}}\right]^N \space |\psi(0)\rangle
    \intertext{Inserting closure and writing in the plane-wave basis: }
     \left| \psi (x, t) \right\rangle &= \int d x_0 \space \langle x_t | e^{\frac{-iĤt}{\hbar}} | x_0 \rangle \space \langle x_0 | \psi(0) \rangle
\end{align*}

Propagation for a single timestep is then:

$$\psi(x,\frac{t_{i+1}}{N}) =
\overbrace{\vphantom{\int \frac{dp}{\sqrt{2\hbar}} e^{\frac{-iV(x)p^{2}}{2 \hbar N}}} e^{\frac{-iV(x)t}{2\hbar N}} }^\textrm{P.E. Propagator} \ \cdot \
\overbrace{\vphantom{\int \frac{dp}{\sqrt{2\hbar}} e^{\frac{-iV(x)p^{2}}{2 \hbar N}}} { \int \frac{dp}{\sqrt{2\pi\hbar}}  \   e^{\frac{-ipx}{\hbar}}} }^\textrm{Inverse Fourier Transform} \ \cdot
\overbrace{\vphantom{\int \frac{dp}{\sqrt{2\hbar}} e^{\frac{-iV(x)p^{2}}{2 \hbar N}}} \ e^{\frac{-ip^2t}{2m\hbar N}}}^\textrm{K.E. Propagator} \ \cdot \
\overbrace{\vphantom{\int \frac{dp}{\sqrt{2\hbar}} e^{\frac{-iV(x)p^{2}}{2 \hbar N}}} {\int \frac{dx}{\sqrt{2\pi\hbar}} \  e^{\frac{ipx}{\hbar}} }}^\textrm{Fourier Transform} \ \cdot \
\overbrace{\vphantom{\int \frac{dp}{\sqrt{2\hbar}} e^{\frac{-iV(x)p^{2}}{2 \hbar N}}}\ e^{\frac{-iV(x)t}{2\hbar N}} }^\textrm{P.E. Propagator} \ \cdot \ \psi(x,t_{i})$$

The Fourier and inverse Fourier transforms are used to convert between the position and momentum basis. To translate the formula above, the algorithm will consist of 5 steps per iteration:

1. Apply a half step of the potential energy propagator to the initial state.
2. Fourier transform into the momentum basis.
3. Apply a full step of the kinetic energy propagator on the momentum basis.
4. Inverse Fourier transform back into the coordinate basis.
5. Apply the second half step of the potential energy propagator.

This is can be done in QFlux with the following code:

```python
# Instantiate our Closed-Systems Dynamics Class
qho_dyn_obj = QFlux_CS(n_basis=128, xo=1.0, po=0.0, mass=1.0, omega=1.0)
# Define our coordinate x and p operators
qho_dyn_obj.set_coordinate_operators()
# Initialize the ladder operators
qho_dyn_obj.initialize_operators()
# Define the default initial state (note that custom initialization is also supported)
qho_dyn_obj.set_initial_state()
# Define some parameters for the time evolution
total_time = 20.0
N_steps = 400
qho_dyn_obj.set_propagation_time(total_time, N_steps)
# Set the Potential/Hamiltonian for our object, in this case using the pre-defined 'harmonic' oscillator potential
qho_dyn_obj.set_hamiltonian(potential_type='harmonic')
# Propagate using the QuTiP sesolve method
qho_dyn_obj.propagate_qt()
# Propagate with SOFT and QuTiP
qho_dyn_obj.propagate_SOFT()
```



## Source Code

::: qflux.closed_systems
    handler: python
    options:
      show_root_heading: true
      show_source: true
      heading_level: 2
      members:
        - DynamicsCS
        - QubitDynamicsCS
        - SpinDynamicsS
        - SpinDynamicsH

