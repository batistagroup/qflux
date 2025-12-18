import numpy as np
import matplotlib.pyplot as plt
from .utils import execute
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit_aer import Aer
from .spin_propagators import get_time_evolution_operator
import numpy.typing as npt
import os 


class SpinDynamicsS:
    """
    A class to simulate the dynamics of a quantum system using a statevector approach.

    Attributes:
        num_qubits (int): Number of qubits in the quantum system.
        evolution_timestep (float): Time step for the evolution.
        trotter_steps (int): Number of Trotter steps for the simulation.
        hamiltonian_coefficients (list): Hamiltonian coefficients for the system.
        initial_state (str or list): Initial state of the system, represented as a binary string or list.
        time_evo_op (QuantumCircuit): Time evolution operator for the system.
        psin0 (ndarray): Initial state vector.
        psin_list (list): List of state vectors during the simulation.
        correlation_list (list): List of correlations calculated during the simulation.
        dpi (int): DPI setting for plot output.
    """

    def __init__(
        self,
        num_qubits,
        evolution_timestep,
        trotter_steps,
        hamiltonian_coefficients,
    ):
        """
        Initialize the QuantumDynamicsClassicalSimulation class.

        Args:
            num_qubits (int): Number of qubits in the quantum system.
            evolution_timestep (float): Time step for the evolution.
            trotter_steps (int): Number of Trotter steps for the simulation.
            hamiltonian_coefficients (list): Hamiltonian coefficients for the system.
        """
        self.num_qubits = num_qubits
        self.evolution_timestep = evolution_timestep
        self.trotter_steps = trotter_steps
        self.hamiltonian_coefficients = hamiltonian_coefficients
        self.time_evo_op = get_time_evolution_operator(
            num_qubits=self.num_qubits,
            tau=self.evolution_timestep,
            trotter_steps=self.trotter_steps,
            coeff=self.hamiltonian_coefficients,
        )
        self.initial_state = ''
        self.psin_list = []
        self.correlation_list = []

        # Plot settings
        self.dpi = 300
        plt.rcParams["axes.linewidth"] = 1.5
        plt.rcParams["lines.markersize"] = 11
        plt.rcParams["figure.figsize"] = (6.4, 3.6)

    def prepare_initial_state(self, state_string):
        """
        Prepare the initial state vector from the binary string or list.

        Returns:
            psin (npt.ArrayLike): Flattened initial state vector.
        """
        zero_state = np.array([[1], [0]])
        one_state = np.array([[0], [1]])

        # Map binary string or list to quantum states
        state_vectors = [
            zero_state if str(bit) == "0" else one_state for bit in state_string
        ]

        # Perform Kronecker product to construct the full initial state
        psin = state_vectors[0]
        for state in state_vectors[1:]:
            psin = np.kron(psin, state)
        return psin.flatten()

    def qsolve_statevector(self, psin):
        """
        Perform statevector propagation for the quantum circuit.
            initial_state (str or list): Initial state of the system, represented as a binary string or list.

        Args:
            psin (npt.ArrayLike): Input statevector.

        Returns:
            ndarray (npt.ArrayLike): Final statevector after execution.
        """
        d = int(np.log2(np.size(psin)))
        qre = QuantumRegister(d)
        circ = QuantumCircuit(qre)
        circ.initialize(psin, qre)
        circ.barrier()
        circ.append(self.time_evo_op, qre)
        circ.barrier()

        device = Aer.get_backend("statevector_simulator")
        result = execute(circ, backend=device).result()
        return result.get_statevector()

    def run_dynamics(self, nsteps, state_string):
        """
        Simulate the dynamics of the quantum system over a number of steps.

        Args:
            nsteps (int): Number of time steps for the simulation.
        """
        self.psin0 = self.prepare_initial_state(state_string)
        self.psin_list = [self.psin0]
        for k in range(nsteps):
            print(f"Running dynamics step {k}")
            if k > 0:
                psin = self.qsolve_statevector(self.psin_list[-1])
                self.psin_list.pop()
                self.psin_list.append(psin)

            self.correlation_list.append(np.vdot(self.psin_list[-1], self.psin0))

    def save_results(self, filename_prefix):
        """
        Save the simulation results to files.

        Args:
            filename_prefix (str): Prefix for the output filenames.
        """
        time = np.arange(
            0,
            self.evolution_timestep * len(self.correlation_list),
            self.evolution_timestep,
        )
        np.save(f"{filename_prefix}_time", time)
        sa_observable = np.abs(self.correlation_list)
        np.save(f"{filename_prefix}_SA_obs", sa_observable)

    def plot_results(self, filename_prefix):
        """
        Plot the simulation results and save the plots as files.

        Args:
            filename_prefix (str): Prefix for the output filenames.
        """
        time = np.arange(
            0,
            self.evolution_timestep * len(self.correlation_list),
            self.evolution_timestep,
        )
        sa_observable = np.abs(self.correlation_list)

        plt.plot(time, sa_observable, "-o")
        plt.xlabel("Time")
        plt.ylabel(r"$\left|\langle \psi | \psi (t)  \rangle \right|$")
        plt.xlim((min(time), max(time)))
        plt.yscale("log")
        #plt.legend()
        plt.tight_layout()
        plt.savefig(f"{filename_prefix}.pdf", format="pdf", dpi=self.dpi)
        plt.savefig(f"{filename_prefix}.png", format="png", dpi=self.dpi)
        plt.show()


class SpinDynamicsH:
    """
    A class to simulate quantum dynamics using Hadamard tests and time evolution operators.

    This class performs quantum circuit construction, execution, and data processing to calculate
    survival amplitudes and probabilities for a spin chain system.
    """
    def __init__(self, num_qubits, evolution_timestep, trotter_steps, hamiltonian_coefficients):
        """
        Initializes the QuantumSimulation class with system parameters and sets up the simulator.

        Parameters:
            num_qubits (int): Number of qubits in the system.
            hamiltonian_coefficients (list): Coefficients for the Hamiltonian terms.
            evolution_timestep (float): Timestep for the time evolution operator.
        """
        self.num_qubits = num_qubits
        self.hamiltonian_coefficients = hamiltonian_coefficients
        self.evolution_timestep = evolution_timestep
        self.trotter_steps = trotter_steps,
        self.simulator = Aer.get_backend('qasm_simulator')
        self.real_amp_list = []
        self.imag_amp_list = []

        self.time_evo_op = get_time_evolution_operator(
            num_qubits=num_qubits,
            tau=evolution_timestep,
            trotter_steps=trotter_steps,
            coeff=hamiltonian_coefficients
        )
        self.controlled_time_evo_op = self.time_evo_op.control()

    def get_hadamard_test(self, initial_state, control_repeats, imag_expectation=False):
        """
        Constructs the Hadamard test circuit to evaluate the real or imaginary components of the operator.

        Parameters:
            initial_state (QuantumCircuit): Circuit for initializing the quantum state.
            control_repeats (int): Number of times the controlled operation is applied.
            imag_expectation (bool): Whether to evaluate the imaginary component (default: False).

        Returns:
            qc (QuantumCircuit): The constructed Hadamard test circuit.
        """
        qr = QuantumRegister(self.num_qubits + 1)
        cr = ClassicalRegister(1)
        qc = QuantumCircuit(qr, cr)
        qc.append(initial_state, qr[1:])
        qc.barrier()

        qc.h(0)
        if imag_expectation:
            qc.p(-np.pi / 2, 0)

        for _ in range(control_repeats):
            qc.append(self.controlled_time_evo_op, qr[:])

        qc.h(0)
        qc.barrier()
        qc.measure(0, 0)
        return qc

    def execute_circuit(self, qc, num_shots=100):
        """
        Executes a quantum circuit using the Qiskit simulator and retrieves measurement counts.

        Parameters:
            qc (QuantumCircuit): The quantum circuit to execute.
            num_shots (int): Number of shots for circuit execution (default: 100).

        Returns:
            dict (dict): Measurement counts.
        """
        job = execute(qc, self.simulator, shots=num_shots)
        return job.result().get_counts()

    @staticmethod
    def calculate_spin_correlation(counts):
        """
        Calculates the spin correlation based on measurement counts.

        Parameters:
            counts (dict): Measurement counts.

        Returns:
            results (float): Average spin correlation.
        """
        qubit_to_spin_map = {'0': 1, '1': -1}
        total_counts = sum(counts.values())
        values = [qubit_to_spin_map[k] * v for k, v in counts.items()]
        return sum(values) / total_counts

    @staticmethod
    def initialize_state(num_qubits, state_string):
        """
        Creates a circuit to initialize the quantum state.

        Parameters:
            num_qubits (int): Number of qubits.
            state_string (str): Binary string representing the initial state.

        Returns:
            qc (QuantumCircuit): Circuit for initializing the state.
        """
        state_string = ''.join('1' if char == '0' else '0' for char in state_string)
        qr = QuantumRegister(num_qubits)
        qc = QuantumCircuit(qr)
        qc.initialize(state_string, qr[:])
        return qc

    def run_simulation(self, state_string, total_time, num_shots=100):
        """
        Runs the Hadamard test simulation for the given initial state.

        Parameters:
            state_string (str): Binary string representing the initial state.
            total_time (float): Total simulation time.
            num_shots (int): Number of shots for circuit execution (default: 100).
        """
        init_circuit = self.initialize_state(self.num_qubits, state_string)
        self.time_range = np.arange(0, total_time + self.evolution_timestep,
                                    self.evolution_timestep)

        for idx, _ in enumerate(self.time_range):
            print(f'Running dynamics step {idx}')

            # Real component
            qc_real = self.get_hadamard_test(init_circuit, idx, imag_expectation=False)
            real_counts = self.execute_circuit(qc_real, num_shots)
            real_amp = self.calculate_spin_correlation(real_counts)
            self.real_amp_list.append(real_amp)

            # Imaginary component
            qc_imag = self.get_hadamard_test(init_circuit, idx, imag_expectation=True)
            imag_counts = self.execute_circuit(qc_imag, num_shots)
            imag_amp = self.calculate_spin_correlation(imag_counts)
            self.imag_amp_list.append(imag_amp)

            print(f'Finished step {idx}: Re = {real_amp:.3f}, Im = {imag_amp:.3f}')

    def save_results(self, prefix):
        """
        Saves the real and imaginary amplitudes, survival amplitude, and survival probability to files.

        Parameters:
            prefix (str): Prefix for the output file names.
        """
        real_amp_array = np.array(self.real_amp_list)
        imag_amp_array = np.array(self.imag_amp_list)

        np.savetxt(f'{prefix}_real_amp.csv', real_amp_array, fmt='%.18e', delimiter=';')
        np.savetxt(f'{prefix}_imag_amp.csv', imag_amp_array, fmt='%.18e', delimiter=';')
        np.savetxt(f'{prefix}_abs_correlation.csv', np.abs(real_amp_array + 1j * imag_amp_array), fmt='%.18e', delimiter=';')
        np.savetxt(f'{prefix}_sqrt_sum_squares.csv', np.sqrt(real_amp_array**2 + imag_amp_array**2), fmt='%.18e', delimiter=';')

    def plot_results(self, prefix):
        """
        Plots the survival amplitude and compares it with a reference.

        Parameters:
            prefix (str): Prefix for loading precomputed reference data.
        """
        abs_corr = np.abs(np.array(self.real_amp_list) + 1j * np.array(self.imag_amp_list))
        
        plt.plot(self.time_range, abs_corr, '.', label='Hadamard Test')
        reference_filename = f'data/{self.num_qubits}_spin_chain_SA_obs.npy'
        if os.path.exists(reference_filename):
            ref_sa = np.load(f'data/{self.num_qubits}_spin_chain_SA_obs.npy')
            ref_time = np.load(f'data/{self.num_qubits}_spin_chain_time.npy')
            plt.plot(ref_time, ref_sa, '-', label='Statevector')
        else:
            print('No reference data file found!')

        plt.xlabel('Time')
        plt.ylabel(r"$\left|\langle \psi | \psi (t)  \rangle \right|$")
        plt.tight_layout()
        #plt.legend()
        plt.show()


if __name__ == "__main__":
    classical = True
    num_q = 3
    evolution_timestep = 0.1
    n_trotter_steps = 1
    hamiltonian_coefficients = [[0.75 / 2, 0.75 / 2, 0.0, 0.65]] + [
        [0.5, 0.5, 0.0, 1.0] for _ in range(num_q - 1)
    ]
    initial_state = "011"  # Specify the initial state as a binary string

    if classical:
        csimulation = SpinDynamicsS(
            num_q,
            evolution_timestep,
            n_trotter_steps,
            hamiltonian_coefficients,
        )
        csimulation.run_dynamics(nsteps=250, state_string=initial_state)
        csimulation.save_results(f"{num_q}_spin_chain")
        csimulation.plot_results(f"{num_q}_spin_chain_statevector")

    else:
        qsimulation = SpinDynamicsH(
            num_q,
            evolution_timestep,
            n_trotter_steps,
            hamiltonian_coefficients,
        )
        qsimulation.run_simulation(state_string=initial_state, total_time=25, num_shots=100)
        qsimulation.save_results('hadamard_test')
        qsimulation.plot_results('hadamard_test')
