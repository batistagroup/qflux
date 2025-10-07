from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit.compiler import transpile
from qiskit_ibm_runtime import QiskitRuntimeService, Sampler
from qiskit.quantum_info.operators import Operator
from qiskit.circuit.library import QFT, PauliEvolutionGate
from qiskit_aer import Aer
from qiskit.synthesis import LieTrotter
import qiskit_aer
import numpy as np
import scipy.linalg as spLA
from tqdm.auto import trange
from .classical_methods import DynamicsCS
from .utils import decompose, pauli_strings_2_pauli_sum
import numpy.typing as npt


class QubitDynamicsCS(DynamicsCS):
    """
    Class to extend `DynamicsCS` by adding qubit-based methods for dynamics.

    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.n_qubits        = int(np.log2(self.n_basis))
        self.quantum_circuit = None


    def _create_QSOFT_Circuit(self, psio: npt.ArrayLike=None):
        """
        Function to construct the QSOFT Circuit.

        Args:
            psio (npt.ArrayLike): initial state that we wish to propagate
        """
        tgrid = self.tlist
        time_step = self.dt
        n_qubits = self.n_qubits
        # Qubit-Basis Propagators
        self.prop_PE_qubit = np.diag(np.exp(-1j*self._PE_grid/2*time_step))
        self.prop_KE_qubit = np.diag(np.exp(-1j*self._KE_grid*time_step))

        q_reg = QuantumRegister(n_qubits)
        c_reg = ClassicalRegister(n_qubits)
        qc = QuantumCircuit(q_reg)
        if type(psio) == type(None):
            qc.initialize(self._psio_grid, q_reg[:], normalize=True)
        else:
            qc.initialize(psio, q_reg[:], normalize=True)
        # Define our PE and KE propagators in Qiskit-friendly manner
        PE_cirq_op = Operator(self.prop_PE_qubit)
        KE_cirq_op = Operator(self.prop_KE_qubit)
        qc.append(PE_cirq_op, q_reg)
        qc.append(QFT(self.n_qubits, do_swaps=True, inverse=False), q_reg)
        qc.append(KE_cirq_op, q_reg)
        qc.append(QFT(self.n_qubits, do_swaps=True, inverse=True), q_reg)
        qc.append(PE_cirq_op, q_reg)
        self.quantum_circuit = qc
        return(qc)


    def _execute_circuit(self, QCircuit: QuantumCircuit, backend=None, shots: int = None, real_backend: bool = False):
        """
            Function to replace the now-deprecated Qiskit
            `QuantumCircuit.execute()` method.

            Args:
                QCircuit (qiskit.QuantumCircuit): qiskit.QuantumCircuit object
                backend (qiskit.Backend): qiskit backend instance
                shots (int): the number of shots to use for circuit sampling

            Returns:
                job: an executed quantum circuit job
        """
        if shots:
            n_shots = shots
        else:
            n_shots = 1024 # Use the qiskit default if not specified
        backend_type = type(backend)
        sv_type = qiskit_aer.backends.statevector_simulator.StatevectorSimulator
        if backend_type == sv_type:
            real_backend = False
        else:
            real_backend = True

        if real_backend:
            QCircuit.measure_all()
            qc = transpile(QCircuit, backend=backend)
            sampler = Sampler(backend)
            job = sampler.run([qc], shots=n_shots)
        else:
            # Transpile circuit with statevector backend
            tmp_circuit = transpile(QCircuit, backend)
            # Run the transpiled circuit
            job = backend.run(tmp_circuit, n_shots=shots)
        return(job)


    def propagate_qSOFT(self, backend=None, n_shots: int = 1024):
        """Function to propagate dynamics object with the qubit SOFT method.

            Args:
                backend (qiskit.Backend): qiskit backend object
                n_shots (int): specifies the number of shots to use when
                    executing the circuit

            Example for using the Statevector Simulator backend:
                >>> from qiskit_aer import Aer
                >>> backend = Aer.get_backend('statevector_simulator')
                >>> self.propagate_qSOFT(backend=backend)
        """
        if backend is None:
            print('A valid backend must be provided ')
        backend_type = type(backend)
        sv_type = qiskit_aer.backends.statevector_simulator.StatevectorSimulator
        if backend_type != sv_type:
            self._propagate_qSOFT_real(backend=backend, n_shots=n_shots)
            return
        else:

            psi_in = self.psio_grid
            # Get initial state from qiskit routine
            q_reg = QuantumRegister(self.n_qubits)
            c_reg = ClassicalRegister(self.n_qubits)
            qc = QuantumCircuit(q_reg, c_reg)
            qc.initialize(self.psio_grid, q_reg[:], normalize=True)
            qc_result = self._execute_circuit(qc, backend=backend, shots=n_shots)
            psio_cirq = qc_result.result().get_statevector().data
            psi_in = psio_cirq
            # Now do propagation loop
            qubit_dynamics_results = [psio_cirq]
            for ii in trange(1, len(self.tlist)):
                circuit = self._create_QSOFT_Circuit(psio=psi_in)
                executed_circuit = self._execute_circuit(circuit, backend=backend, shots=n_shots)
                psi_out = executed_circuit.result().get_statevector().data
                qubit_dynamics_results.append(psi_out)
                psi_in = psi_out

            self.dynamics_results_qSOFT = np.asarray(qubit_dynamics_results)
            return


    def get_statevector_from_counts(self, counts, n_shots):
        new_statevector = np.zeros_like(self.psio_grid)

        for key in counts:
            little_endian_int = int(key, 2)
            new_statevector[little_endian_int] = counts[key]/n_shots
        return(new_statevector)


    def _propagate_qSOFT_real(self, backend='statevector_simulator', n_shots=1024):
        """
            Function to propagate dynamics object with the qubit SOFT method.

            Args:
                backend (qiskit.Backend): qiskit backend object
                n_shots (int): specifies the number of shots to use when
                    executing the circuit

            Example for using the Statevector Simulator backend:
                >>> from qiskit_aer import Aer
                >>> backend = Aer.get_backend('statevector_simulator')
                >>> self.propagate_qSOFT(backend=backend)
        """


        psi_in = self.psio_grid
        # Get initial state from qiskit routine
        q_reg = QuantumRegister(self.n_qubits)
        c_reg = ClassicalRegister(self.n_qubits, name='c')
        qc = QuantumCircuit(q_reg, c_reg)
        qc.initialize(self.psio_grid, q_reg[:], normalize=True)
        # Now do propagation loop
        qubit_dynamics_results = []
        for ii in trange(len(self.tlist)):
            circuit = self._create_QSOFT_Circuit(psio=psi_in)
            executed_circuit = self._execute_circuit(circuit, backend=backend, shots=n_shots)
            circuit_result = executed_circuit.result()
            measured_psi = circuit_result[0].data['meas'].get_counts()
            self._last_measurement = measured_psi
            psi_out = self.get_statevector_from_counts(measured_psi, n_shots)
            psi_in = psi_out
            qubit_dynamics_results.append(psi_out)
            psi_in = psi_out

        self.dynamics_results_qSOFT = np.asarray(qubit_dynamics_results)
        return


    def _construct_pauli_gate(self, hamiltonian_matrix=None):
        """
            Function to construct a pauli evolution gate from Hamiltonian
        
            Args:
                hamiltonian_matrix (npt.ArrayLike): array-like matrix representing the hamiltonian of interest
                    If not provided, use the operator representation of the Hamiltonian by default.

        """

        if type(hamiltonian_matrix) == type(None):
            decomposed_H = decompose(self.H_op.full())
        else:
            decomposed_H = decompose(hamiltonian_matrix)
        H_pauli_sum  = pauli_strings_2_pauli_sum(decomposed_H)
        synthesizer  = LieTrotter(reps=2)
        prop_pauli_H = PauliEvolutionGate(operator=H_pauli_sum, time=self.dt, synthesis=synthesizer)
        self.pauli_prop = prop_pauli_H
        self._pauli_hamiltonian = decomposed_H
        return

    
    def _construct_pauli_cirq(self, psio=None):
        q_reg = QuantumRegister(self.n_qubits)
        c_reg = ClassicalRegister(self.n_qubits)
        qc = QuantumCircuit(q_reg, c_reg)
        qc.initialize(psio, q_reg[:], normalize=True)
        qc.append(self.pauli_prop, q_reg)
        self.quantum_circuit = qc
        return(qc)


    def propagate_qmatvec(self, backend=None, n_shots: int = 1024, hamiltonian_matrix=None, initial_state=None):
        r"""
            Function to propagate dynamics object with the qubit matvec method.

            Args:
                backend (qiskit.Backend): qiskit backend object
                n_shots (int): specifies the number of shots to use when
                    executing the circuit
                hamiltonian_matrix (npt.ArrayLike): array-like matrix representing the Hamiltonian
                    Used to construct the propagator:

                    $$ U(t) = e^{- i H t / \hbar} $$ 

                    By default, the operator representation of the hamiltonian `self.H_op` is used.
                initial_state (npt.ArrayLike): array-like vector representing the initial state
            Example for using the Statevector Simulator backend:
                >>> from qiskit_aer import Aer
                >>> backend = Aer.get_backend('statevector_simulator')
                >>> self.propagate_qSOFT(backend=backend)
        """

        # Create the Pauli propagator:
        self._construct_pauli_gate(hamiltonian_matrix=hamiltonian_matrix)

        q_reg = QuantumRegister(self.n_qubits)
        c_reg = ClassicalRegister(self.n_qubits)
        qc = QuantumCircuit(q_reg, c_reg)
        # Initialize State
        if type(initial_state) == type(None):
            qc.initialize(self.psio_op.full().flatten(), q_reg[:], normalize=True)
        else:
            qc.initialize(initial_state, q_reg[:], normalize=True)
        qc_result = self._execute_circuit(qc, backend=backend, shots=n_shots)
        psio_cirq = qc_result.result().get_statevector().data
        psi_in = psio_cirq
        new_qubit_dynamics_result = [psio_cirq]
        for ii in trange(1, len(self.tlist)):
            circuit = self._construct_pauli_cirq(psio=psi_in)
            executed_circuit = self._execute_circuit(circuit, backend=backend, shots=n_shots)
            psi_out = executed_circuit.result().get_statevector().data
            new_qubit_dynamics_result.append(psi_out)
            psi_in = psi_out
        self.dynamics_results_qubit = np.asarray(new_qubit_dynamics_result)
        return

