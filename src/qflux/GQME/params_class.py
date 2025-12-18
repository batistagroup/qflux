"""
Example Parameters class as illustrated in Part_I_JCE.ipynb.

"""

import numpy as np

class Params:
    def __init__(self):
        # ==== Spin-Boson Model parameters ====
        self.GAMMA_DA = 1          # diabatic coupling
        self.EPSILON  = 1
        self.BETA     = 5          # inverse temperature beta = 1 / (k_B * T)
        self.XI       = 0.1
        self.OMEGA_C  = 2

        # Spin-up and spin-down states
        self.spin_up   = np.array([1.0, 0.0], dtype=np.float64)
        self.spin_down = np.array([0.0, 1.0], dtype=np.float64)

        # ==== General constants for simulation ====
        self.TIME_STEPS = 500                    # number of time steps
        self.au2ps      = 0.00002418884254       # as -> a.u. conversion
        self.timeau     = 12.409275
        self.DT         = 20 * self.au2ps * self.timeau  # time step in au

        self.FINAL_TIME = self.TIME_STEPS * self.DT
        self.DOF_E      = 2                      # number of electronic states
        self.DOF_E_SQ   = self.DOF_E * self.DOF_E

        # ==== Simulation parameters for TT-TFD ====
        self.DOF_N      = 50     # number of nuclear DOF
        self.OMEGA_MAX  = 10

        # TT constants
        self.eps         = 1e-12         # tt approx error
        self.dim         = self.DOF_N    # number of coords
        self.occ         = 10            # max occupation number
        self.MAX_TT_RANK = 10

        # ==== Simulation parameters for GQME ====
        self.MEM_TIME          = self.DT * self.TIME_STEPS
        self.HBAR              = 1
        self.MAX_ITERS         = 30
        self.CONVERGENCE_PARAM = 10.0**(-10.0)

        # ==== Parameter string for output files ====
        self.PARAM_STR  = "_Spin-Boson_Ohmic_TT-TFD_b%sG%s_e%s_" % (
            self.BETA, self.GAMMA_DA, self.EPSILON
        )
        self.PARAM_STR += "xi%swc%s_wmax%s_dofn%s" % (
            self.XI, self.OMEGA_C, self.OMEGA_MAX, self.DOF_N
        )

        # ==== Pauli matrices ====
        self.X = np.array([[0, 1], [1, 0]], dtype=np.complex128)
        self.Y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
        self.Z = np.array([[1, 0], [0, -1]], dtype=np.complex128)
        self.I = np.eye(2, dtype=np.complex128)
