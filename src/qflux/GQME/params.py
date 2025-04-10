import numpy as np

## Spin-Boson Model parameters 
GAMMA_DA = 1 # diabatic coupling
EPSILON = 1
BETA = 5 # inverse finite temperature beta = 1 / (k_B * T)
XI = 0.1
OMEGA_C = 2

print("SPIN-BOSON Model parameter")
print("        epsilon =", EPSILON)
print("             xi =", XI)
print("        omega_c =", OMEGA_C)

# Spin-up and spin-down states
spin_up = np.array([1.0, 0.0], dtype=np.float64)
spin_down = np.array([0.0, 1.0], dtype=np.float64)

## General Constants for simulation
TIME_STEPS = 500          # number of time steps
au2ps = 0.00002418884254   # Conversion of attoseconds to atomic units
timeau = 12.409275
DT = 20 * au2ps * timeau  # time step in au

FINAL_TIME = TIME_STEPS * DT # final time
DOF_E = 2 # number of electronic states
DOF_E_SQ = DOF_E * DOF_E

##Simulation Parameter for TT-TFD
# TFD parameter: for Discretized nuclear DOFs
DOF_N = 50 # number of nuclear DOF
OMEGA_MAX = 10

# TT constants
eps = 1e-12           # tt approx error
dim = DOF_N         # number of coords
occ = 10                 # maximum occupation number; low for harmonic systems
MAX_TT_RANK = 10

print("      omega_max =", OMEGA_MAX)
print("     time steps =", TIME_STEPS)
print("             DT =", DT)
print("     final time =", FINAL_TIME)
print("          DOF_E =", DOF_E)
print("          DOF_N =", DOF_N)

##Simulation Parameter for GQME
MEM_TIME = DT * TIME_STEPS
HBAR = 1
MAX_ITERS = 30
CONVERGENCE_PARAM = 10.**(-10.)

##Parameter for output files
PARAM_STR = "_Spin-Boson_Ohmic_TT-TFD_b%sG%s_e%s_"%(BETA, GAMMA_DA, EPSILON)
PARAM_STR += "xi%swc%s_wmax%s_dofn%s"%(XI, OMEGA_C, OMEGA_MAX, DOF_N)

# Pauli matrices
X = np.array([[0, 1], [1, 0]], dtype=np.complex128)
Y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
Z = np.array([[1, 0], [0, -1]], dtype=np.complex128)
I = np.eye(2, dtype=np.complex128)

