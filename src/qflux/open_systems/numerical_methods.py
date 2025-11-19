import numpy as np
import scipy.linalg as LA
import scipy.fft as sfft
from qutip import mesolve, Qobj
from typing import List, Tuple, Callable, Optional, Any

from . import params as pa
from . import trans_basis as tb


class DynamicsOS:
    """Class for open-system dynamics (Lindblad equation).

    This class provides methods to simulate open-system dynamics described by the Lindblad equation.

    Attributes:
        Nsys (int): System Hilbert Space Dimension.
        Hsys (np.ndarray): Hamiltonian of the system (shape (N, N)).
        rho0 (np.ndarray): Initial density matrix (shape (N, N)).
        c_ops (List[np.ndarray]): List of collapse operators (each of shape (N, N)).
    """

    def __init__(
        self,
        Nsys: int,
        Hsys: np.ndarray,
        rho0: np.ndarray,
        c_ops: Optional[List[np.ndarray]] = None
    ) -> None:
        """
        Initialize the DynamicsOS instance.

        Args:
            Nsys (int): System Hilbert Space Dimension.
            Hsys (np.ndarray): Hamiltonian of the system.
            rho0 (np.ndarray): Initial density matrix.
            c_ops (Optional[List[np.ndarray]]): List of collapse operators. Defaults to an empty list.
        """
        if c_ops is None:
            c_ops = []
        self.Nsys: int = Nsys
        self.Hsys: np.ndarray = Hsys
        self.rho0: np.ndarray = rho0
        self.c_ops: List[np.ndarray] = c_ops

    def Gt_matrix_expo(self, time_arr: List[float], Is_show_step: bool = False) -> List[np.ndarray]:
        """
        Compute the propagator of the Lindblad equation using the matrix exponential.

        The propagator is computed by exponentiating the Liouvillian operator defined by the
        system Hamiltonian and collapse operators.

        Args:
            time_arr (List[float]): Array of time values for the simulation.
            Is_show_step (bool, optional): If True, prints the current simulation step. Defaults to False.

        Returns:
            List[np.ndarray]: List of propagators corresponding to each time in `time_arr`.
        """
        ident_h: np.ndarray = np.eye(self.Nsys, dtype=np.complex128)

        # Build the A matrix for time-derivation of the vectorized density matrix.
        Amat: np.ndarray = -1j * (np.kron(self.Hsys, ident_h) - np.kron(ident_h, self.Hsys.T))
        for i in range(len(self.c_ops)):
            Amat += 0.5 * (
                2.0 * np.kron(self.c_ops[i], self.c_ops[i].conj())
                - np.kron(ident_h, (self.c_ops[i].T @ self.c_ops[i].conj()))
                - np.kron((self.c_ops[i].T.conj() @ self.c_ops[i]), ident_h)
            )

        G_prop: List[np.ndarray] = []
        for i, t in enumerate(time_arr):
            if Is_show_step:
                print("step", i, "time", t)
            Gt: np.ndarray = LA.expm(Amat * t)
            G_prop.append(Gt)
        return G_prop

    def propagate_matrix_exp(
        self,
        time_arr: List[float],
        observable: np.ndarray,
        Is_store_state: bool = False,
        Is_show_step: bool = False,
        Is_Gt: bool = False,
    ) -> Any:
        """
        Solve the Lindblad equation using matrix exponential.

        This method computes the propagator, evolves the initial density matrix, and calculates
        the expectation value of the observable over time. Optionally, it stores the evolved density matrices.

        Args:
            time_arr (List[float]): Time array for dynamic simulation.
            observable (np.ndarray): Observable matrix for which the expectation value is computed.
            Is_store_state (bool, optional): If True, stores the density matrices at each time step.
                                             Defaults to False.
            Is_show_step (bool, optional): If True, prints the current simulation step. Defaults to False.
            Is_Gt (bool, optional): If True, includes the propagators in the result. Defaults to False.

        Returns:
            Result: An object with the following attributes:
                - expect (List[float]): List of expectation values over time.
                - density_matrix (List[np.ndarray], optional): List of density matrices (if `Is_store_state` is True).
                - Gprop (List[np.ndarray], optional): List of propagators (if `Is_Gt` is True).
        """

        class Result:
            """Class for storing propagation results."""
            def __init__(self, store_state: bool, include_Gt: bool) -> None:
                self.expect: List[float] = []
                if store_state:
                    self.density_matrix: List[np.ndarray] = []
                if include_Gt:
                    self.Gprop: Optional[List[np.ndarray]] = None

        result = Result(Is_store_state, Is_Gt)

        # Compute the propagator of the Lindblad equation.
        G_prop: List[np.ndarray] = self.Gt_matrix_expo(time_arr, Is_show_step)
        if Is_Gt:
            result.Gprop = G_prop

        # Initialize the vectorized density matrix.
        vec_rho0: np.ndarray = self.rho0.reshape(self.Nsys**2)

        for i, _ in enumerate(time_arr):
            vec_rhot: np.ndarray = G_prop[i] @ vec_rho0
            # Reshape back to density matrix form.
            rhot: np.ndarray = vec_rhot.reshape(self.Nsys, self.Nsys)

            if Is_store_state:
                result.density_matrix.append(rhot)
            result.expect.append(np.trace(rhot @ observable).real)

        return result

    def propagate_qt(self, time_arr: List[float], observable: Any, **kwargs: Any) -> List[float]:
        """
        Propagate the system using QuTiP's `mesolve` function.

        This method solves the Lindblad master equation using QuTiP's `mesolve` to compute the expectation
        values of the observable over time.

        Args:
            time_arr (List[float]): Time array for dynamic simulation.
            observable (Any): Observable operator(s) for which the expectation value is computed.
                              Can be a single operator or a list of operators.
            **kwargs: Additional keyword arguments to pass to `mesolve`.

        Returns:
            List[float]: List of expectation values of the observable over time.
        """
        c_ops: List[Qobj] = [Qobj(c_op) for c_op in self.c_ops]

        if isinstance(observable, list):
            obs = [Qobj(op) for op in observable]
        else:
            obs = Qobj(observable)

        result = mesolve(Qobj(self.Hsys), Qobj(self.rho0), time_arr, c_ops=c_ops, e_ops=obs, **kwargs)
        return result.expect


class DVR_grid:
    """Class for Discrete Variable Representation (DVR) grid methods.

    This class handles grid-based representations for systems where the potential is expressed on grid points.

    Attributes:
        Ngrid (int): Number of grid points.
        xmin (float): Minimum value of the grid.
        xmax (float): Maximum value of the grid.
        mass (float): Mass of the particle.
        xgrid (np.ndarray): Array of grid points in position space.
        dx (float): Spacing between grid points.
        dk (float): Spacing in momentum space.
        kgrid (np.ndarray): Array of grid points in momentum space.
        ak2 (np.ndarray): Kinetic energy array in momentum space.
        hamk (np.ndarray): Kinetic Hamiltonian matrix in position space.
        potential (Optional[np.ndarray]): Potential energy array on the grid.
    """

    def __init__(self, xmin: float, xmax: float, Ngrid: int, mass: float) -> None:
        """
        Initialize the DVR_grid instance.

        Args:
            xmin (float): Minimum x-value.
            xmax (float): Maximum x-value.
            Ngrid (int): Number of grid points.
            mass (float): Mass of the particle.
        """
        self.Ngrid: int = Ngrid
        self.xmin: float = xmin
        self.xmax: float = xmax
        self.mass: float = mass

        # Set up the position grid.
        self.xgrid: np.ndarray = np.array([])
        self._set_xgrid()
        self.dx: float = self.xgrid[1] - self.xgrid[0]

        # Set up the momentum grid.
        self.dk: float = 2.0 * np.pi / (self.Ngrid * self.dx)
        self.kgrid: np.ndarray = np.array([])
        self.ak2: np.ndarray = np.array([])  # Kinetic energy array.
        self.hamk: np.ndarray = np.array([])  # Kinetic Hamiltonian matrix.
        self._set_kinet_ham()

        # Potential energy array (to be set later).
        self.potential: Optional[np.ndarray] = None

    def _set_xgrid(self) -> None:
        """
        Set up the position space grid.

        Initializes the `xgrid` attribute using a linear space between `xmin` and `xmax`.
        """
        self.xgrid = np.linspace(self.xmin, self.xmax, self.Ngrid)

    def set_potential(self, func_pot: Callable[[float], float]) -> None:
        """
        Set up the potential energy array on the grid.

        Args:
            func_pot (Callable[[float], float]): Function that returns the potential value at a given x.
        """
        self.potential = np.zeros_like(self.xgrid)
        for i in range(self.Ngrid):
            self.potential[i] = func_pot(self.xgrid[i])

    def _set_kinet_ham(self) -> None:
        """
        Set up the kinetic Hamiltonian matrix in position space.

        This method computes the momentum grid and the corresponding kinetic energy values,
        and constructs the kinetic Hamiltonian matrix in position space using a Fourier transform.
        """
        self.kgrid = np.zeros(self.Ngrid, dtype=np.float64)
        self.ak2 = np.zeros(self.Ngrid, dtype=np.float64)

        coef_k: float = pa.hbar**2 / (2.0 * self.mass)

        for i in range(self.Ngrid):
            if i < self.Ngrid // 2:
                self.kgrid[i] = i * self.dk
            else:
                self.kgrid[i] = -(self.Ngrid - i) * self.dk
            self.ak2[i] = coef_k * self.kgrid[i]**2

        akx0: np.ndarray = sfft.ifft(self.ak2)
        self.hamk = np.zeros((self.Ngrid, self.Ngrid), dtype=np.complex128)

        for i in range(self.Ngrid):
            for j in range(self.Ngrid):
                if i < j:
                    self.hamk[i, j] = akx0[i - j].conj()
                else:
                    self.hamk[i, j] = akx0[i - j]

    def get_eig_state(self, Nstate: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get the eigenstates for the potential in x-space.

        Args:
            Nstate (int): Number of eigenstates to output.

        Returns:
            Tuple[np.ndarray, np.ndarray]: A tuple containing:
                - Eigenvalues (np.ndarray) for the first `Nstate` states.
                - Eigenvectors (np.ndarray) for the first `Nstate` states, normalized by sqrt(dx).
        """
        Mata: np.ndarray = self.hamk.copy()
        for i in range(self.Ngrid):
            Mata[i, i] += self.potential[i]

        val, arr = LA.eigh(Mata)
        return val[:Nstate], arr[:, :Nstate] / np.sqrt(self.dx)

    def x2k_wave(self, psi: np.ndarray) -> np.ndarray:
        """
        Transform the wavefunction from position space to momentum space.

        Args:
            psi (np.ndarray): Wavefunction in position space.

        Returns:
            np.ndarray: Wavefunction in momentum space.
        """
        return tb.x2k_wave(self.dx, psi)

    def k2x_wave(self, psik: np.ndarray) -> np.ndarray:
        """
        Transform the wavefunction from momentum space to position space.

        Args:
            psik (np.ndarray): Wavefunction in momentum space.

        Returns:
            np.ndarray: Wavefunction in position space.
        """
        return tb.k2x_wave(self.dx, psik)
