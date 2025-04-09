# **Simulation of a Spin-Chain System**

In this example, we simulate a spin-chain system comprised of three spins. We can do things first using a classical Statevector simulation, then we can run things in the Quantum Circuit formalism. 

## Statevector Simulation

```python
from qflux.closed_systems.spin_propagators import * 
from qflux.closed_systems.hamiltonians import * 
from qflux.closed_systems.spin_dynamics_oo import * 

num_q = 3
evolution_timestep = 0.1
n_trotter_steps = 1
hamiltonian_coefficients = [[0.75 / 2, 0.75 / 2, 0.0, 0.65]] + [[0.5, 0.5, 0.0, 1.0]
                            for _ in range(num_q - 1)]
initial_state = "011"  # Specify the initial state as a binary string

csimulation = SpinDynamicsS(
                            num_q,
                            evolution_timestep,
                            n_trotter_steps,
                            hamiltonian_coefficients
                            )
csimulation.run_dynamics(nsteps=250, state_string=initial_state)
csimulation.save_results(f"{num_q}_spin_chain")
csimulation.plot_results(f"{num_q}_spin_chain_statevector")

```

## Quantum Circuit Simulation with Hadamard Test

```python
num_q = 3
evolution_timestep = 0.1
n_trotter_steps = 1
hamiltonian_coefficients = [[0.75 / 2, 0.75 / 2, 0.0, 0.65]] + [[0.5, 0.5, 0.0, 1.0]
                            for _ in range(num_q - 1)]
initial_state = "011"  # Specify the initial state as a binary string

qsimulation = SpinDynamicsH(
                            num_q,
                            evolution_timestep,
                            n_trotter_steps,
                            hamiltonian_coefficients
                            )
qsimulation.run_simulation(state_string=initial_state, total_time=25, num_shots=100)
qsimulation.save_results('hadamard_test')
qsimulation.plot_results('hadamard_test')
```

