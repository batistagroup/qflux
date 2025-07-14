import numpy as np
import scipy.linalg as LA
from typing import List, Tuple, Optional, Dict, Any, Union

from qiskit import transpile
from qiskit_aer import AerSimulator
#from qiskit.primitives import Estimator
from qiskit_ibm_runtime import EstimatorV2 as Estimator
from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit import QuantumCircuit

from . import trans_basis as tb
from . import dilation_circuit as dc
from .numerical_methods import DynamicsOS


def expand(Gmat_org: np.ndarray, Norg: int, Nexpand: int) -> np.ndarray:
    """
    Expand the propagator in the vectorized density matrix representation.

    This function expands the original propagator matrix to match the dimensions
    of the vectorized density matrix representation (e.g., match the 2^N dimension).

    Args:
        Gmat_org (np.ndarray): Original propagator matrix.
        Norg (int): Original system dimension.
        Nexpand (int): Expanded system dimension.

    Returns:
        np.ndarray: The expanded propagator matrix with shape (Nexpand**2, Nexpand**2).
    """
    Gnew = np.zeros((Nexpand**2, Nexpand**2), dtype=np.complex128)
    for i in range(Norg):
        for j in range(Norg):
            for k in range(Norg):
                for l in range(Norg):
                    Gnew[i * Nexpand + j, k * Nexpand + l] = Gmat_org[i * Norg + j, k * Norg + l]
    return Gnew


def gen_Kraus_list(Gmat: np.ndarray, N: int, tol: float = 1e-5) -> List[np.ndarray]:
    """
    Generate the Kraus operators from the propagator with a given tolerance.

    This function computes the Kraus operator representation by constructing the Choi
    matrix from the propagator matrix and performing an eigenvalue decomposition.

    Args:
        Gmat (np.ndarray): Matrix of the propagator with shape (N^2, N^2).
        N (int): The system Hilbert space dimension.
        tol (float, optional): Tolerance threshold for the Kraus operator representation.
            Operators with eigenvalues below tol are ignored. Defaults to 1e-5.

    Returns:
        List[np.ndarray]: A list of Kraus operators.
    """
    # Define the Choi matrix from the propagator matrix
    C_mat = np.zeros(Gmat.shape, dtype=np.complex128)
    for i in range(N):
        for j in range(N):
            C_matij = np.zeros(Gmat.shape, dtype=np.complex128)
            for k in range(N):
                for l in range(N):
                    C_matij[i * N + k, l * N + j] = Gmat[j * N + k, l * N + i]
            C_mat += C_matij

    Kraus: List[np.ndarray] = []
    eigenvalues, eigenvectors = LA.eigh(C_mat)
    for idx in range(len(eigenvalues)):
        if eigenvalues[idx] > tol:
            Mi = np.sqrt(eigenvalues[idx]) * eigenvectors[:, idx].reshape(N, N)
            Kraus.append(Mi.conj().T)
    return Kraus


class QubitDynamicsOS(DynamicsOS):
    """
    Class for simulating quantum dynamics using either vectorized density matrix or Kraus operator representations.

    This class provides methods to initialize state vectors, construct quantum circuits,
    and perform quantum simulations using Qiskit backends.
    """

    def __init__(self, rep: str = 'Density', **kwargs: Any) -> None:
        """
        Initialize a QubitDynamicsOS instance.

        Depending on the representation, either "Density" or "Kraus", the number of qubits is computed.
        Additional keyword arguments are passed to the base DynamicsOS class.

        Args:
            rep (str, optional): Representation type, either 'Density' or 'Kraus'. Defaults to 'Density'.
            **kwargs: Additional keyword arguments for the DynamicsOS initializer.
        """
        super().__init__(**kwargs)

        if rep == 'Density':
            # Vectorized density matrix representation
            self.rep: str = 'Density'
            self.Nqb: int = int(np.log2(self.Nsys**2))
        elif rep == 'Kraus':
            # Kraus operator representation
            self.rep = 'Kraus'
            self.Nqb = int(np.log2(self.Nsys))

        # The counting qubits bit string and observable matrix are initialized to None.
        self.count_str: Optional[List[str]] = None
        self.observable: Optional[np.ndarray] = None

        # Default dilation method for quantum simulation.
        self.dilation_method: str = 'Sz-Nagy'

    def set_dilation_method(self, method: str) -> None:
        """
        Set the dilation method for quantum simulation.

        Args:
            method (str): The dilation method, e.g., 'Sz-Nagy', 'SVD', or 'SVD-Walsh'.
        """
        self.dilation_method = method

    def set_count_str(self, count_str: List[str]) -> None:
        """
        Set the counting bit string for measurement.

        Args:
            count_str (List[str]): The counting bit string.
        """
        self.count_str = count_str

    def set_observable(self, observable: np.ndarray) -> None:
        """
        Set the observable for the quantum simulation.

        Args:
            observable (np.ndarray): The observable matrix.
        """
        self.observable = observable

    def init_statevec_vecdens(self) -> Tuple[np.ndarray, float]:
        """
        Initialize the state vector from the initial density operator using vectorized representation.

        The initial density matrix is reshaped into a vector and normalized.

        Returns:
            Tuple[np.ndarray, float]: A tuple containing the normalized state vector and the norm
            of the original vectorized density matrix.
        """
        vec_rho0 = self.rho0.reshape(self.Nsys**2)
        norm0 = LA.norm(vec_rho0, 2)
        statevec = vec_rho0 / norm0
        return statevec, norm0

    def init_statevec_Kraus(self, tol: float = 1e-6) -> Tuple[List[np.ndarray], List[float]]:
        """
        Initialize state vectors from the initial density operator using the Kraus operator representation.

        The density matrix is decomposed using eigenvalue decomposition, and eigenstates with eigenvalues
        below the specified tolerance are ignored.

        Args:
            tol (float, optional): Tolerance for ignoring eigenstates with small eigenvalues. Defaults to 1e-6.

        Returns:
            Tuple[List[np.ndarray], List[float]]: A tuple containing a list of state vectors and a list of
            corresponding probabilities.
        """
        eigenvalues, eigenvectors = LA.eigh(self.rho0)
        statevec: List[np.ndarray] = []
        prob: List[float] = []
        for i in range(len(eigenvalues) - 1, -1, -1):
            if abs(eigenvalues[i]) < tol:
                break
            prob.append(eigenvalues[i])
            statevec.append(eigenvectors[:, i])
        return statevec, prob

    def _get_qiskit_observable(self, Isdilate: bool = False, tol: float = 5e-3) -> SparsePauliOp:
        """
        Prepare and return the Qiskit observable operator.

        Converts the observable matrix to its Pauli representation and returns a SparsePauliOp.

        Args:
            Isdilate (bool, optional): Flag indicating whether to use the dilated observable.
                Defaults to False.
            tol (float, optional): Tolerance for the Pauli decomposition. Defaults to 5e-3.

        Returns:
            SparsePauliOp: The Qiskit representation of the observable.
        """
        if self.observable is None:
            print('Error: observable is None')

        if Isdilate:
            num_qubits = self.Nqb + 1
            Obs_mat = np.zeros((2 * self.Nsys, 2 * self.Nsys), dtype=np.complex128)
            Obs_mat[:self.Nsys, :self.Nsys] = self.observable[:self.Nsys, :self.Nsys]
        else:
            num_qubits = self.Nqb
            Obs_mat = self.observable

        Obs_paulis_dic = tb.ham_to_pauli(Obs_mat, num_qubits, tol=tol)

        # Prepare the Qiskit observable from the Pauli strings of the observable matrix.
        data: List[str] = []
        coef: List[float] = []
        for key in Obs_paulis_dic:
            data.append(key)
            coef.append(Obs_paulis_dic[key])
        obs_q = SparsePauliOp(data, coef)
        return obs_q

    def qc_simulation_kraus(
        self,
        time_arr: List[float],
        Kraus: Optional[Dict[int, List[np.ndarray]]] = None,
        backend: Any = AerSimulator(), 
        Gprop: Optional[List[np.ndarray]] = None,
        tolk: float = 1e-5,
        tolo: float = 1e-5,
        Is_store_circuit = False,
        **kwargs: Any
        ) -> np.ndarray:
        """
        Perform quantum simulation using the Kraus operator representation.

        This method simulates the quantum system dynamics over a series of time steps using a Kraus operator-based approach.
        It constructs quantum circuits for each Kraus operator and initial state, runs the simulation using Qiskit's Estimator,
        and accumulates the measurement results.

        Args:
            time_arr (List[float]): List of time steps for simulation.
            
            Kraus (Optional[Dict[int, List[np.ndarray]]], optional): Dictionary mapping time step index to a list of Kraus operators.
                If None, Kraus operators are generated from the propagator. Defaults to None.
            backend (Any, optional): Quantum simulation backend. Defaults to AerSimulator().
            Gprop (Optional[List[np.ndarray]], optional): Propagator matrix (or list of matrices) for simulation.
                If None, it will be calculated. Defaults to None.
            tolk (float, optional): Tolerance for generating Kraus operators. Defaults to 1e-5.
            tolo (float, optional): Tolerance for observable decomposition. Defaults to 1e-5.
            Is_store_circuit (bool, optional): If True, store the generated quantum circuits at each time step. Defaults to False.
            **kwargs: Additional keyword arguments for propagator calculation.

        Returns:
            Dict[str, Any]: A dictionary containing quantum simulation results, including:
                - 'data' (np.ndarray): Array of shape (nsteps,), containing the accumulated observable expectation values at each time step.
                - 'circuits' (List[List[QuantumCircuit]], optional): If `Is_store_circuit` is True, this field contains a list of lists of quantum circuits.
                  Each outer list corresponds to a time step, and each inner list contains all circuits used for that step (across different Kraus operators and initial states).
        """
        nsteps = len(time_arr)

        # Generate Kraus operators if not provided.
        if Kraus is None:
            Kraus = {}
            if Gprop is None:
                print('Calculating the propagator')
                Gprop = self.Gt_matrix_expo(time_arr, **kwargs)
            print('Generating the Kraus operators')
            for i in range(nsteps):
                print('At step', i, 'of', nsteps)
                Kraus[i] = gen_Kraus_list(Gprop[i], self.Nsys, tol=tolk)
        print('Kraus operator generation complete')

        # Perform Qiskit simulation using the Estimator.
        estimator = Estimator(mode=backend)

        statevec, prob = self.init_statevec_Kraus()
        n_inistate = len(statevec)
        print('Number of initial states in the density matrix:', n_inistate)
        print('Probabilities:', prob)

        # Obtain the Qiskit observable.
        obs_q = self._get_qiskit_observable(Isdilate=True, tol=tolo)

        print('Starting quantum simulation')
        
        result_simu = {}
        result_simu['data'] = np.zeros(nsteps, dtype=np.float64)
        if(Is_store_circuit):  result_simu['circuits'] = [ [] for _ in range(nsteps) ]

        for i in range(nsteps):
            print('Simulation step', i, 'of', nsteps)
            current_kraus_list = Kraus[i]
            print('Number of Kraus operators:', len(current_kraus_list))
            for kraus_op in current_kraus_list:
                for istate in range(n_inistate):
                    qc = self._create_circuit(kraus_op, statevec[istate], Isscale=False)
                    result = estimator.run([(qc, obs_q)]).result()
                    
                    if(Is_store_circuit):  result_simu['circuits'][i].append(qc)
                    result_simu['data'][i] += result[0].data.evs * prob[istate]

        return result_simu

    def qc_simulation_vecdens(
        self,
        time_arr: List[float],
        shots: int = 1024,
        backend: Any = AerSimulator(),
        Gprop: Optional[List[np.ndarray]] = None,
        **kwargs: Any
    ) -> np.ndarray:
        """
        Perform quantum simulation using the vectorized density matrix representation.

        This method simulates the quantum system dynamics over a series of time steps by constructing circuits
        based on the vectorized density matrix representation, performing measurements, and processing the results.

        Args:
            time_arr (List[float]): List of time steps for simulation.
            shots (int, optional): Number of measurement shots. Defaults to 1024.
            backend (Any, optional): Quantum simulation backend. Defaults to AerSimulator().
            Gprop (Optional[List[np.ndarray]], optional): Propagator matrix (or list of matrices) for simulation.
                If None, it will be calculated. Defaults to None.
            **kwargs: Additional keyword arguments for propagator calculation.

        Returns:
            Dict[str, Any]: A dictionary containing the quantum simulation results, including:
                - 'data' (np.ndarray): Array of shape (nsteps, n_bitstr), storing the processed probability amplitudes 
                (i.e., normalized square root of the measured probabilities, scaled by norm factors) at each time step.
                - 'circuits' (List[QuantumCircuit]): List of generated quantum circuits for each time step.
                - 'norm' (List[float]): List of norm correction factors used in post-processing at each time step.
        """
        if Gprop is None:
            Gprop = self.Gt_matrix_expo(time_arr, **kwargs)

        nsteps = len(time_arr)

        if self.count_str is None:
            print("Error: count_str is not assigned")

        n_bitstr = len(self.count_str)
        statevec, norm0 = self.init_statevec_vecdens()
        
        result = {}
        result['data'] = np.zeros((nsteps, n_bitstr), dtype=np.float64)
        result['circuits'] = []
        result['norm'] = []
        
        for i in range(nsteps):
            if i % 100 == 0:
                print('Quantum simulation step', i)
            Gt = Gprop[i]
            
            #create the circuit
            circuit, norm = self._create_circuit(Gt, statevec, Isscale=True)
            circuit.measure(range(self.Nqb + 1), range(self.Nqb + 1))
            if self.dilation_method == 'SVD-Walsh':
                circuit = transpile(circuit, backend)
            
            #store the circuits and norm to the result
            result['circuits'].append(circuit)
            result['norm'].append(norm)
            
            counts = backend.run(circuit, shots=shots).result().get_counts()
            for j in range(n_bitstr):
                bitstr = self.count_str[j]
                if bitstr in counts:
                    result['data'][i, j] = np.sqrt(counts[bitstr] / shots) * norm * norm0
                else:
                    print('At time', i, 'with shots =', shots, "no counts for", bitstr)

        return result

    def _create_circuit(
        self,
        array: np.ndarray,
        statevec: Union[np.ndarray, List[np.ndarray]],
        Isscale: bool = True
    ) -> QuantumCircuit:
        """
        Construct and return the quantum circuit.

        This method wraps the call to the dilation circuit construction function.

        Args:
            array (np.ndarray): Array used for circuit construction (e.g., propagator or Kraus operator).
            statevec (Union[np.ndarray, List[np.ndarray]]): State vector(s) to be used in the circuit.
            Isscale (bool, optional): Flag indicating whether scaling should be applied. Defaults to True.

        Returns:
            QuantumCircuit: The constructed quantum circuit.
        """
        return dc.construct_circuit(self.Nqb, array, statevec, method=self.dilation_method, Isscale=Isscale)
