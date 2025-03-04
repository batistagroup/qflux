"""
Define some Parameter and Constants
"""

import numpy as np

# Pauli matrices
X = np.array([[0, 1], [1, 0]], dtype=np.complex128)
Y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
Z = np.array([[1, 0], [0, -1]], dtype=np.complex128)
I = np.eye(2, dtype=np.complex128)


# sigma+ and sigma-
sigmap = 0.5*(X+1j*Y)
sigmam = 0.5*(X-1j*Y)

# Spin-up and spin-down states
spin_up = np.array([1.0, 0.0], dtype=np.float64)
spin_down = np.array([0.0, 1.0], dtype=np.float64)


#some constant to convert the units 
au2fs    = 2.418884254E-2
au2cm    = 219474.63068
au2joule = 4.35974381E-18
bolz     = 1.3806503E-23
au2ev    = 27.2114
hbar     = 1.0