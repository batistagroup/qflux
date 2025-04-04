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

### Anatomy of a Closed System Quantum Dynamics Simulation

To run a dynamics simulation, we need to define some key quantities. This section will walk you through the process of defining these things within `qflux`.

As a reminder, the task at hand is compute the time evolution of a wavefunction according to the Schrodinger equation: 

$$ \left| \psi (t) \right\rangle = e^{- \frac{i}{\hbar} H t} \left| \psi_{0} \right\rangle $$ 

To do this, we must: 

- Define the initial state $\psi_{0}$. 
- Define the Hamiltonian describing the system of interest.
- Define a propgation time-step $t$ and the number of time steps $n$ for which we want to compute the evolved wavefunction. 

We will now look at how each of these steps can be done with qflux. 

#### Definition of the Initial State

The initial state is our wavefunction $\psi_{0}$. In order to define this abstract object on a computer, we must define a finite space in which it exists. The number of discrete points in this space is controlled by the `n_basis` parameter that is passed to the `Dynamics_CS` class upon instantiation. Note that if you do not define this argument, the default value of 128 is used. Given 128 grid points, we can begin defining operators. To compute the dynamics in the so-called "Fock basis", we define the ladder operators $\hat{a}$, $\hat{x}$, and $\hat{p}$. This is done by calling the `.intialize_operator()` method. To compute the dynamics in the position/coordinate basis, we must define a range of position-values that define the x-grid of our space. This can be done with the `.set_coordinate_operators(x_min=-7., x_max=7)` method, which will define an array of `n_basis` points, ranging from `x_min` to `x_max`. Now that we've defined the space in which our wavefunction can exist, we can finally define the wavefunction. 

When instantiating a dynamics object with the `Dynamics_CS` class, there are some other important arguments that are taken into account: 

- `xo`: The initial displacement in the position-coordinate. 
- `po`: The initial displacement in the momentum-coordinate. 
- `mass`: The mass of the particle/system of interest.

Note that these should all be defined in atomic units. 

To define the initial state in our default way, you can simply use the `.set_initial_state()` method. This takes the optional argument of `wfn_omega` defining the frequency/width of the intiial state, which takes the default value of 1.0 au. The default initial state in the coordinate basis is defined as a Gaussian coherent state: 

$$ \psi_{0} = \left( \frac{m \omega}{\pi \hbar} \right)^{1/4} e^{- \frac{m \omega}{2 \hbar} \left( x - x_{0} \right)^{2} + \frac{i}{\hbar} p_{0} x} $$ 

The default initial state in the ladder/Fock basis is defined as the coherent state with amplitude $\alpha = (x_{0} + i p_{0})/ \sqrt{2}$, defined in the Fock basis as: 

$$ \left| \alpha \right\rangle = \sum_{n=0}^{n_{basis}} \frac{\alpha^{n}}{\sqrt{n!}} e^{- \frac{1}{2} \left| \alpha \right|^{2} \left| n \right\rangle $$ 

where $n$ is a state in the Fock basis. 

Note that `qflux` also provides functionality for custom state initialization, in which a user-defined function can be provided. 

Custom initialization in the coordinate basis is done with the `.custom_grid_state_initialization()` method. This can be used if you wish to initialize with a different state. To exemplify the usage of this functionality, we can define some function that takes arguments necessary to define a state: 

```python
import numpy as np

def custom_gaussian(xvals, xo, po, omega, mass, hbar):
    normalization_factor = (mass * omega / (np.pi * hbar))**(0.25)
    exp_func = np.exp( - mass * omega / (2 * hbar) * (xvals - xo)**2 + 1j/hbar * po * xvals)
    return(normalization_factor * exp_func)
```

Then, we can set-up a dynamics object that is ready to initialize a state: 

```python 
HO_dyn_obj = Dynamics_CS(n_basis=128, xo=1.0, po=0.0, mass=1.0, omega=0.2)
HO_dyn_obj.set_coordinate_operators()
HO_dyn_obj.initialize_operators()
```

And define some arguments to provide to the function: 

```python
func_args = {'xvals': HO_dyn_obj.x_grid,
             'xo'   : HO_dyn_obj.xo,
             'po'   : HO_dyn_obj.po,
             'omega': 1.0,
             'mass': HO_dyn_obj.mass,
             'hbar': 1.0}
```

And we can now define the state by calling: 

```python
HO_dyn_obj.custom_grid_state_initialization(custom_gaussian, **func_args)
```


Similarly, to define a custom initial state to be used in the Fock/ladder basis, we follow a similar pattern. We first define some custom function that will return a `qutip.Qobj`: 

```python
def custom_coherent_state(N_basis):
    '''
    Function to define a squeezed coherent state by squeezing, then displacing the vacuum state.
    '''
    squeezed_coh_state = qt.displace(N_basis, 2) * qt.squeeze(N_basis, 1.0) * qt.basis(N_basis, 0)
    return(squeezed_coh_state)
```

And define the keyword arguments as a dictionary:

```python
qt_func_args = {'N_basis': HO_dyn_obj.n_basis}
```

Then we can perform the custom initialization by calling the `.custom_ladder_state_initialization()` method:

```python
HO_dyn_obj.custom_ladder_state_initialization(custom_coherent_state, **qt_func_args)
```

For the custom Fock/ladder basis initialization, the custom function must return a `qutip.Qobj`. 

#### Defining the Hamiltonian

The next step for running our dynamics simulation is to define the Hamiltonian, which should describe the system of interest. For the coordinate basis, we assume that the Hamiltonian takes the form of: 

$$ H = V(x) + \frac{p^{2}}{2 m} $$ 

where $V(x)$ describes the potential energy of our system. `qflux` provides some example systems out of the box, which we will now demonstrate how to use.

To use the built-in potential energy functions, all one must do is use the `.set_hamiltonian()` method. This method has the optional keyword argument `potential_type`, which can be used to choose one of the two currently implemented potentials:
- Harmonic Oscillator Potential
- Arbitrary Quartic Potential


The harmonic oscilator potential is implemented in the grid-basis as: 

$$ V(x) = \frac{1}{2} m \omega^{2} x^{2} $$ 

and in the ladder basis as: 

$$ H = \hbar \omega \left( \hat{a}^{\dagger} \hat{a} + \frac{1}{2} \right) $$ 

The frequency ($\omega$) and mass ($m$) can be controlled when instantiating the dynamics object with the `mass` and `omega` keyword arguments. 


The quartic potential is implemented as: 

$$ V(x) = a_{0} + a_{1} \frac{x}{x_{0}} + a_{2} \frac{x^{2}}{x_{0}^{2}} + a_{3} \frac{x^{3}}{x_{0}^{3}} + a_{4} \frac{x^{4}}{x_{0}^{4}} $$ 

To use a custom quartic potential, the user should provide a dictionary of keyword arguments that define the coefficients ($a_{0}, a_{1}, a_{2}, a_{3}, a_{4}$) and the scaling factor $x_{0}$: 

```python
coeffs_dict = {'a0': 1, 'a1': 1, 'a2': 1, 'a3': 1, 'a4': 1, 'x0': 1}
dyn_obj.set_hamiltonian(potential_type='quartic', **coeffs_dict)
```

For the Fock/ladder basis, the $x$ in the previous equation is replaced with an operator $\hat{a}$ defined in terms of the creation and annihilation operators as $\hat{x} = \frac{1}{\sqrt{2}} \left( \hat{a}^{\dagger} + \hat{a} \right)$ and in the kinetic energy term $\frac{p^{2}}{2m}$, $p$ is replaced with $\hat{p} =  \frac{i}{\sqrt{2}} \left( \hat{a}^{dagger} - \hat{a} \right)$.

`qflux` also supports arbitrary customization of the potential energy function by use of the `.set_H_grid_with_custom_potential()` and `.set_H_op_with_custom_potential()` methods. These methods expect a function and a dictionary (of keyword arguments for that function) as arguments. This is illustrated in the following example: 

Suppose you have some arbitrary Morse-like potential of the form:

$$ V_{Morse} = De ( 1 - e^{- a (x-x_{eq}))^{2} $$

We can define a python function to to construct this potential:

```python
def morse_potential(x_eq=None, mass=None, omega=None, xval=None):
    De = 8
    xe = 0
    k = mass*omega**2
    a = np.sqrt(k/(2*De))
    y = De * ((1 - np.exp(-a*(xval-x_eq)))**2)
    return(y)
```

And define a dictionary with parameters to define a specific potential:

```python
morse_args = {'x_eq': -1.0, 'mass': 1.0, 'omega': 1.5, 'xval': dyn_obj.x_grid}
```

And then construct a Hamiltonian with this custom function for our dynamics object by calling: 

```python
dyn_obj.set_H_grid_with_custom_potential(morse_potential, **morse_args)
```

Similarly, we can do this in a `qutip.Qobj`-compatible format: 

```python
def morse_potential_op(x_eq=None, mass=None, omega=None, xval=None):
    De = 8
    xe = 0
    k = mass*omega**2
    a = np.sqrt(k/(2*De))
    exponential_f = (-a * (xval - x_eq)).expm()
    y = De * ((1 - exponential_f)**2)
    return(y)
```

And define our dictionary of custom parameters: 

```python
morse_op_args = {'x_eq': -1.0, 'mass': 1.0, 'omega': 1.5, 'xval': dyn_obj.x_op}
```

And finally we can set our Hamiltonian for the Fock/ladder basis by calling: 

```python
dyn_obj.set_H_op_with_custom_potential(morse_potential_op, **morse_op_args)
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

