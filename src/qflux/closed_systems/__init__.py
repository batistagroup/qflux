"""

qflux v.0.1.0

__version__ = 0.1.0

Description: 

classical_methods.py --> DynamicsCS : Closed-System dynamics 
qubit_methods.py     --> QubitDynamicsCS : Qubit-based dynamics
spin_dynamics_oo.py  --> SpinDynamicsS, SpinDynamicsH : Statevector and hadamard-test implementations (respectively) for spin-chain systems

"""

from .classical_methods import DynamicsCS
from .qubit_methods import QubitDynamicsCS
from .spin_dynamics_oo import SpinDynamicsS, SpinDynamicsH
from .direct_method import hamiltonian_simulation
