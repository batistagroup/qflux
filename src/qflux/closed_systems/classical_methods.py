# Class for classical propagation methods
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import axes
from tqdm.auto import trange
import qutip as qt
from typing import Callable
from .utils import *


class DynamicsCS:
    """
        Class for closed-system dynamics. **All input parameters must be in
        atomic units to ensure consistency. Please be sure to convert your
        parameters to atomic units prior to instantiation.**

    """

    def __init__(self, n_basis: int = 128, xo: float = 1.0, po: float = 0.0, mass: float = 1.0,
                 omega: float = 1.0) -> None:
        """
        Args:
            n_basis (int): Number of states to include in the chosen representation. If basis
                = 'ladder', this is the Fock cutoff and defines the number of states
                used for representing the ladder operators. If basis = 'coordinate',
                this defines the number of points for the position and momenta.

            xo (float, optional): Defines the displacement of the initial state in the position
                coordinate. Default is 1.0 Bohr.

            po (float, optional): Defines the displacement of the initial state in the position
                coordinate. Default is 1.0 au.

            mass (float, optional): Defines the mass of the particle/system of interest.
                Default is 1.0 au.

            omega (float, optional): Frequency of harmonic oscillator.
                Default is 1.0 au.

        """
        #--------- Required Attributes Populated During Execution ----------#
        self.n_basis                   = n_basis
        self.xo                        = xo
        self.po                        = po
        self.mass                      = mass
        self.hbar                      = 1.0
        self.omega                     = omega
        #--------- Below are Attributes Populated During Execution ---------#
        self.total_time                = 0.
        self.n_tsteps                  = 0.
        self._KE_op                    = None
        self._PE_op                    = None
        self.H_op                      = None
        self.prop_KE_op                = None
        self.prop_PE_op                = None
        self.prop_H_op                 = None
        # Grid operators
        self._KE_grid                  = None
        self._PE_grid                  = None
        self.H_grid                    = None
        self.PE_prop_grid              = None
        self.KE_prop_grid              = None


    def _get_xgrid(self, x_min: float, x_max: float) -> None:
        """
        Populate the `self.x_grid` and `self.dx` attributes. This function
        generates an array of `self.n_basis` evenly spaced values between
        `x_min` and `x_max`.

        Args:
            x_min (float): Minimum value of x-coordinates
            x_max (float): Maximum value of x-coordinates

        Returns:
            self.dx (float): Spacing between points in the x-coordinate grid.
            self.xgrid (array_like): Array of grid points from x_min to x_max with spacing of dx
        """
        dx = (x_max - x_min) / self.n_basis
        x_grid = np.arange(-self.n_basis / 2, self.n_basis / 2) * dx
        self.dx = dx
        self.x_grid = x_grid
        return


    def _get_pgrid(self, x_min: float, x_max: float, reorder: bool = True) -> None:
        """
        Populate the `self.p_grid` and `self.dp` attributes. This function
        generates an array of `self.n_basis` evenly spaced values.

        Args:
            x_min (float): Minimum value of x-coordinates
            x_max (float): Maximum value of x-coordinates
            reorder (bool): Boolean flag to determine whether points should be reordered to be
                compatible with the FFT routine or not.

        Returns:
            self.dp (float): Spacing between points in the p-coordinate grid.
            self.pgrid (array_like): Array of momentum grid points
        """
        dp = 2 * np.pi / (x_max - x_min)
        pmin = -dp * self.n_basis / 2
        pmax = dp * self.n_basis / 2
        plus_pgrid = np.linspace(0, pmax, self.n_basis // 2 + 1)
        minus_pgrid = - np.flip(np.copy(plus_pgrid))
        if reorder:
            p_grid = np.concatenate((plus_pgrid[:-1], minus_pgrid[:-1]))
        else:
            p_grid = np.concatenate((minus_pgrid, plus_pgrid))
        self.p_grid = p_grid
        self.dp = dp
        return


    def set_coordinate_operators(self, x_min: float = -7., x_max: float = 7., reorder_p: bool = True) -> None:
        """
        Populate the `self.x_grid`, `self.p_grid`, `self.dx`, and `self.dp`
        attributes. This functions generates an array of `self.n_basis`
        evenly spaced values.

        Args:
            x_min : float
                Minimum value of x-coordinates
            x_max : float
                Maximum value of x-coordinates
            reorder_p : bool
                Boolean flag to determine whether momentum values should be
                reordered to be compatible with the FFT routine or not.

        Returns:
            self.dx : float
                Spacing between points in the x-coordinate grid.
            self.xgrid : array_like
                Array of x-values
            self.dp : float
                Spacing between points in the p-coordinate grid.
            self.pgrid : array_like
                Array of p-values
        """
        self._get_xgrid(x_min, x_max)
        self._get_pgrid(x_min, x_max, reorder=reorder_p)
        return


    def initialize_operators(self):
        """
            Function to initialize core operators in the chosen basis.

        """

        self.a_op = qt.destroy(self.n_basis)
        self.x_op = qt.position(self.n_basis)
        self.p_op = qt.momentum(self.n_basis)
        return


    def _set_hamiltonian_grid(self, potential_type: str = 'harmonic', **kwargs):
        if potential_type == 'harmonic':

            # Set attributes for the coordinate basis
            self._PE_grid = self.mass * self.omega ** 2 * self.x_grid ** 2 / 2.
            self._KE_grid = self.p_grid ** 2 / (2. * self.mass)
            self.H_grid = self._PE_grid + self._KE_grid
        elif potential_type == 'quartic':
            if kwargs:
                a0 = 0.0
                a1 = 0.0
                a2 = 0.0
                a3 = 0.0
                a4 = 0.0
                if 'a0' in kwargs:
                    a0 = kwargs['a0']
                if 'a1' in kwargs:
                    a1 = kwargs['a1']
                if 'a2' in kwargs:
                    a2 = kwargs['a2']
                if 'a3' in kwargs:
                    a3 = kwargs['a3']
                if 'a4' in kwargs:
                    a4 = kwargs['a4']
                if 'x0' in kwargs:
                    x0 = kwargs['x0']
                # Assume that all inputs have the proper atomic units:
                cf = 1.0
                xi = self.x_op

            else:
                # Define relevant parameters
                cf = convert_eV_to_au(1.)
                x0 = 1.9592
                a0 = 0.0
                a1 = 0.429
                a2 = -1.126
                a3 = -0.143
                a4 = 0.563
                # Do calculation for ladder basis
                xi = self.x_grid / x0
            self._PE_grid = cf * (a0 + a1 * xi + a2 * xi ** 2 + a3 * xi ** 3 + a4 * xi ** 4)
            self._KE_grid = self.p_grid ** 2 / (2. * self.mass)
            self.H_grid = self._PE_grid + self._KE_grid
        return


    def _set_hamiltonian_qt(self, potential_type: str = 'harmonic', **kwargs):
        if potential_type == 'harmonic':
            # Set attributes for the ladder basis
            self.H_op = self.hbar * self.omega * (self.a_op.dag() * self.a_op + 0.5)
            self._KE_op = self.p_op ** 2 / (2. * self.mass)
            self._PE_op = self.mass * self.omega ** 2 * self.x_op ** 2 / 0.5
            self.H_xp_op = self._PE_op + self._KE_op
        elif potential_type == 'quartic':
            if kwargs:
                a0 = 0.0
                a1 = 0.0
                a2 = 0.0
                a3 = 0.0
                a4 = 0.0
                if 'a0' in kwargs:
                    a0 = kwargs['a0']
                if 'a1' in kwargs:
                    a1 = kwargs['a1']
                if 'a2' in kwargs:
                    a2 = kwargs['a2']
                if 'a3' in kwargs:
                    a3 = kwargs['a3']
                if 'a4' in kwargs:
                    a4 = kwargs['a4']
                if 'x0' in kwargs:
                    x0 = kwargs['x0']
                # Assume that all inputs have the proper atomic units:
                cf = 1.0
                xi = self.x_op

            else:
                # Define relevant parameters
                cf = convert_eV_to_au(1.)
                x0 = 1.9592
                a0 = 0.0
                a1 = 0.429
                a2 = -1.126
                a3 = -0.143
                a4 = 0.563
                # Do calculation for ladder basis
                xi = self.x_op / x0
            self.x0 = x0
            self._PE_op = cf * (a0 + a1 * xi + a2 * xi ** 2 + a3 * xi ** 3 + a4 * xi ** 4)
            self._KE_op = self.p_op ** 2 / (2. * self.mass)
            self.H_op = self._PE_op + self._KE_op
            return

    def set_hamiltonian(self, potential_type: str = 'harmonic', **kwargs):
        """
        Function to define Hamiltonian.

        Args:
            potential_type : str
                String defining the type of potential energy surface.
                Available options are: ('harmonic', 'quartic', ...)

                Note: You can manually define your potential energy using the functions:
                    - set_H_grid_with_custom_potential
                    - set_H_op_with_custom_potential

        """

        if potential_type == 'harmonic':
            self._set_hamiltonian_grid(potential_type=potential_type, **kwargs)
            self._set_hamiltonian_qt(potential_type=potential_type, **kwargs)
        elif potential_type == 'quartic':
            self._set_hamiltonian_grid(potential_type=potential_type, **kwargs)
            self._set_hamiltonian_qt(potential_type=potential_type, **kwargs)
        else:
            print('Error, this potential type has not yet been implemented!')
            print('Set your parameters with the custom functions!')
        return


    def set_H_grid_with_custom_potential(self, custom_function: Callable, **kwargs):
        """
        Function to allow for user-defined potential defined by custom_function. Must be a function of qutip operators.

        Args:
            custom_function (Callable): Function that defines the custom potential
                energy. Must return an array

        """
        potential = custom_function(**kwargs)
        self._PE_grid = potential
        self._KE_grid = self.p_grid ** 2 / (2. * self.mass)
        self.H_grid = self._PE_grid + self._KE_grid
        return


    def set_H_op_with_custom_potential(self, custom_function: Callable, **kwargs):
        """
        Function to allow for user-defined potential defined by custom_function. Must be a function of qutip operators.

        Args:
            custom_function (Callable): Function that defines the potential
                energy in terms of qutip QObj operators. Must return a qutip.Qobj
        """
        potential = custom_function(**kwargs)
        self._PE_op = potential
        self._KE_op = self.p_op ** 2 / (2. * self.mass)
        self.H_op = self._PE_op + self._KE_op
        return


    def set_initial_state(self, wfn_omega: float = 1.0):
        """
        Function to define the initial state. By default, a coherent state is
        used as the initial state defined in the basis chosen upon instantiation

        Args:
            wfn_omega (float, optional): Defines the frequency/width of the initial state.
                Default is 1.0 au.
        """

        alpha_val = (self.xo + 1j * self.po) / np.sqrt(2)
        psio = qt.coherent(self.n_basis, alpha=alpha_val)
        # Now populate the initial state in the grid basis
        normalization = (self.mass * wfn_omega / np.pi / self.hbar) ** (0.25)
        exponential = np.exp(-1 * (self.mass * wfn_omega / self.hbar / 2) *
                             ((self.x_grid - self.xo) ** 2)
                             + 1j * self.po * self.x_grid / self.hbar
                             )

        coherent_state = normalization * exponential
        # Set the attributes
        self.psio_grid = coherent_state
        self.psio_op = psio
        return



    def custom_grid_state_initialization(self, function_name: Callable, **kwargs):
        """
        Function to allow for customized grid state initialization.

        Args:
            function_name (Callable): name of user-defined function that returns
                the initial state. Must return an array
        """

        self.psio_grid = function_name(**kwargs)
        return

    def custom_ladder_state_initialization(self, function_name: Callable, **kwargs):
        """
        Function to allow for customized ladder state initialization.

        Args:
            function_name (Callable): name of user-defined function that returns
                the initial state. Must return a qutip.Qobj.
        """

        self.psio_op = function_name(**kwargs)
        return

    def set_propagation_time(self, total_time: float, n_tsteps: int):
        """
        Function to define the propagation time, an array of times from
        t=0 to total_time, with n_tsteps equally-spaced steps.

        Args:
        total_time : float
            The total time for which we wish to compute the dynamics.
        n_tsteps : int
            The number of equally-spaced time steps used to compute the dynamics

        Returns:
        self.tlist : array-like

        """

        self.tlist = np.linspace(0., total_time, n_tsteps+1)
        self.dt = self.tlist[1] - self.tlist[0]
        self.n_tsteps = n_tsteps
        return


    def propagate_qt(self, solver_options : dict = None):
        """
        Function used to propagate with qutip.

        Args:
            solver_options (dict): A dictionary of arguments to pass to the qutip.sesolve function

        Returns:
            dynamics_results (array-like): array containing the propagated state

        """

        options = {'nsteps': len(self.tlist),
                    'progress_bar': True}

        if solver_options:
            for key in solver_options:
                options[key] = solver_options[key]

        results = qt.sesolve(self.H_op, self.psio_op, self.tlist,
                             options=options)

        self.dynamics_results_op = results.states
        return


    def propagate_SOFT(self):
        """
        Function used to propagate with the 2nd-Order Trotter Expansion.

        $$
        e^{- \\frac{i}{\\hbar} H t} \\approx e^{- \\frac{i}{\\hbar} V t/2} e^{- \\frac{i}{\\hbar} T t} e^{- \\frac{i}{\\hbar} V t/2} + \\mathcal{O}^{3}
        $$

        Returns:
            dynamics_results_grid (array-like): array containing the propagated state
                                                shape (n_tsteps x self.n_basis)

        """
        self.tau = self.tlist[1] - self.tlist[0]
        PE_prop = np.exp(-1.0j * self._PE_grid / 2 * self.tau / self.hbar)
        KE_prop = np.exp(-1.0j * self._KE_grid * self.tau / self.hbar)

        self.PE_prop_grid = PE_prop
        self.KE_prop_grid = KE_prop

        propagated_states = [self.psio_grid]
        psi_t = self.psio_grid
        for ii in range(1, len(self.tlist)):
            psi_t_position_grid = PE_prop * psi_t
            psi_t_momentum_grid = KE_prop * np.fft.fft(psi_t_position_grid, norm="ortho")
            psi_t = PE_prop * np.fft.ifft(psi_t_momentum_grid, norm="ortho")
            propagated_states.append(psi_t)

        self.dynamics_results_grid = np.asarray(propagated_states)
        return
