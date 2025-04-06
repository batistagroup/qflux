# Anatomy of a Closed System Quantum Dynamics Simulation

To run a dynamics simulation, we need to define some key quantities. This section will walk you through the process of defining these things within `qflux`.

As a reminder, the task at hand is compute the time evolution of a wavefunction according to the Schrodinger equation: 

$$ \left| \psi (t) \right\rangle = e^{- \frac{i}{\hbar} H t} \left| \psi_{0} \right\rangle $$ 

To do this, we must: 

- Define the initial state $\psi_{0}$. 
- Define the Hamiltonian describing the system of interest.
- Define a propgation time-step $t$ and the number of time steps $n$ for which we want to compute the evolved wavefunction. 

We will now look at how each of these steps can be done with qflux. 

## Definition of the Initial State

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
dyn_obj = Dynamics_CS(n_basis=128, xo=1.0, po=0.0, mass=1.0, omega=0.2)
dyn_obj.set_coordinate_operators()
dyn_obj.initialize_operators()
```

And define some arguments to provide to the function: 

```python
func_args = {'xvals': dyn_obj.x_grid,
             'xo'   : dyn_obj.xo,
             'po'   : dyn_obj.po,
             'omega': 1.0,
             'mass': dyn_obj.mass,
             'hbar': 1.0}
```

And we can now define the state by calling: 

```python
dyn_obj.custom_grid_state_initialization(custom_gaussian, **func_args)
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
qt_func_args = {'N_basis': dyn_obj.n_basis}
```

Then we can perform the custom initialization by calling the `.custom_ladder_state_initialization()` method:

```python
dyn_obj.custom_ladder_state_initialization(custom_coherent_state, **qt_func_args)
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

$$ V_{Morse} = De ( 1 - e^{- a (x-x_{eq})})^{2} $$

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


