import numpy as np
import numpy.typing as npt
from typing import List, Optional, Tuple

from qiskit import QuantumCircuit
from qiskit_aer.primitives import EstimatorV2 as Estimator
from qiskit.quantum_info import SparsePauliOp
from qiskit_aer.noise import NoiseModel
from qiskit_ibm_runtime.fake_provider import FakeSherbrooke

# To change the ansatz, apply_param and measure_der must both be modified.
def apply_param(
    params: npt.NDArray[np.float64], i: int, qc: QuantumCircuit, N: int
) -> None:
    """Apply parameter i to the quantum circuit currently constructing the ansatz.
        The ansatz must be built in a peicewise manner to allow for hadamard tests
        of the generators of the parameters to later be inserted to measure the A_ij
        and C_i matrix elements.

    Args:
        params (numpy.array): An array containing the values of all the ansatz parameters.
        parameter (int): Index of the parameter being applied.
        qc (QuantumCircuit): The qiskit ansatz quantum circuit currently being constructed.
        N (int): Number of qubits
    """
    qc.rx(params[i], i % N)
    if i % N == N - 1 and i != len(params) - 1:
        for i in range(N - 1):
            qc.cz(i, i + 1)


def measure_der(i: int, qc: QuantumCircuit, N: int) -> None:
    """Append a Hadamard test to the circuit to measure the generator of parameter i in the ansatz.
        The ansatz currently used is simply the two-local ansatz with only rx gates.
        Therefore the generator is only x-gates on the corresponding qubit.

    Args:
        parameter (int): The index of the parameter whose generator will be measured.
        qc (QuantumCircuit): The qiskit quantum circuit which is in the process of assembling the ansatz.
        N (int): Number of qubits
    """
    qc.cx(N, i % N)


def A_Circuit(params: npt.NDArray[np.float64], i: int, j: int, N: int) -> QuantumCircuit:
    """Constructs the qiskit quantum circuits used to measure each element of the A_ij matrix.

    Args:
        params (numpy.array): A numpy array containing the values of each parameter of the ansatz.
        i (int): The index from A_ij.  This also corresponds to the ansatz parameter i being measured.
        j (int): The index from A_ij.  This also corresponds to the ansatz parameter j being measured.
        N (int): The number of qubits.

    Returns:
        QuantumCircuit: The quantum circuit for an element of the A_ij matrix, in the form of a
            hadamard test of the generators of parameters i and j.  The value of the A_ij matrix
            can be found by measuring the ancilla qubit (qubit N) in the Z basis.
    """
    qc = QuantumCircuit(N + 1, 1)
    qc.h(N)
    for parameter in range(len(params)):  # Apply parameterized gates
        if parameter == i:
            qc.x(N)
            measure_der(parameter, qc, N)  # Measure generator for i
            qc.x(N)
        if parameter == j:
            measure_der(parameter, qc, N)  # Measure second generator for j
        apply_param(params, parameter, qc, N)
    qc.h(N)
    return qc


def Measure_A(
    init_circ: QuantumCircuit,
    params: npt.NDArray[np.float64],
    N: int,
    shots: int = 2**10,
    noisy: bool = False,
) -> npt.NDArray[np.float64]:
    """Create the A_ij matrix through measuring quantum circuits corresponding to each element.

    Args:
        init_circ (QuantumCircuit): The qiskit circuit representing the initial state of the system.
        params (numpy.array): A numpy array which contains the values of each parameter of the ansatz.
        N (int): The number of qubits.
        shots (int, optional): The number of shots used to estimate each element of the A_ij matrix. Defaults to 2**10.
        noisy (bool, optional): A boolean used to turn on and off the Fake-Sherbrooke qiskit noisy backend. Defaults to False.

    Returns:
        numpy.array: The A_ij matrix
    """
    A = [[0.0 for i in range(len(params))] for j in range(len(params))]
    for i in range(len(params)):
        for j in range(len(params) - i):
            qc = QuantumCircuit(N + 1, 1)
            ansatz = A_Circuit(params, i, i + j, N)
            qc = qc.compose(init_circ, [k for k in range(N)])
            qc = qc.compose(ansatz, [k for k in range(N + 1)])

            observable = SparsePauliOp.from_list([("Z" + "I" * N, 1.0)])
            if noisy:
                device_backend = FakeSherbrooke()
                coupling_map = device_backend.coupling_map
                noise_model = NoiseModel.from_backend(device_backend)
                basis_gates = noise_model.basis_gates
                estimator = Estimator(options={
                                    "backend_options":{"noise_model": noise_model},
                                    "run_options":{"shots": shots}}
                )
            else:
                estimator = Estimator(options={"run_options":{"shots": shots}})
            result = estimator.run([(qc, observable)]).result()
            A[i][i + j] = result[0].data.evs
    return np.array(A)


def C_Circuit(
    params: npt.NDArray[np.float64],
    i: int,
    pauli_string: str,
    N: int,
    evolution_type: str = "real",
) -> QuantumCircuit:

    """Create the qiskit quantum circuits to measure each element of the C_i vector.

    Args:
        params (numpy.array): A numpy array which contains the values of each parameter of the ansatz.
        i (int): The index of the C_i vector being measured. This also corresponds
            to the index i of the parameter whose generator will be measured
        pauli_string (str): A string containing a description of the pauli operator of the Hamiltonian which will be measured.
        N (int): The number of qubits.
        evolution_type (str, optional): This determines if the evolution will be real-time or imaginary-time
            through the addition of an extra gate. Defaults to "real".

    Returns:
        QuantumCircuit: The quantum circuit for an element of the C_i matrix, in the form of a
            hadamard test of the generators of parameter i.  The value of the C_i matrix
            can be found by measuring the ancilla qubit (qubit N) in the Z basis.
    """
    qc = QuantumCircuit(N + 1, 1)
    qc.h(N)
    if evolution_type == "imaginary":
        qc.s(N)  # To get only imaginary component
    else:
        qc.z(N)
    for parameter in range(len(params)):  # Apply parameterized gates
        if parameter == i:
            qc.x(N)
            measure_der(parameter, qc, N)  # Measure generators
            qc.x(N)
        apply_param(params, parameter, qc, N)
    pauli_measure(qc, pauli_string)
    qc.h(N)
    return qc


def Measure_C(
    init_circ: QuantumCircuit,
    params: npt.NDArray[np.float64],
    H: SparsePauliOp,
    N: int,
    shots: int = 2**10,
    evolution_type: str = "real",
    noisy: bool = False,
) -> npt.NDArray[np.float64]:
    """Create the C_i vector through measuring quantum circuits corresponding to each element.

    Args:
        init_circ (QuantumCircuit): A qiskit circuit constructing the initial state of the system.
        params (numpy.array): A numpy array containing the values of the parameters of the ansatz.
        H (SparsePauliOp): The Hamiltonian.
        N (int): The number of qubits.
        shots (int, optional): The number of shots to be used to measure each element of the C_i vector. Defaults to 2**10.
        evolution_type (str, optional): This determines if the evolution will be real-time or imaginary-time
            through the addition of an extra gate. Defaults to "real".
        noisy (bool, optional): A boolean used to turn on and off the Fake-Sherbrooke qiskit noisy backend. Defaults to False.

    Returns:
        numpy.array: The C_i vector.
    """
    C = [0.0 for i in range(len(params))]
    for i in range(len(params)):
        for pauli_string in range(len(H.paulis)):
            qc = QuantumCircuit(N + 1, 1)
            ansatz = C_Circuit(
                params, i, H.paulis[pauli_string], N, evolution_type=evolution_type
            )
            qc = qc.compose(init_circ, [k for k in range(N)])
            qc = qc.compose(ansatz, [k for k in range(N + 1)])
            observable = SparsePauliOp.from_list([("Z" + "I" * N, 1.0)])
            if noisy:
                device_backend = FakeSherbrooke()
                coupling_map = device_backend.coupling_map
                noise_model = NoiseModel.from_backend(device_backend)
                basis_gates = noise_model.basis_gates
                estimator = Estimator(options={
                                        "backend_options":{"noise_model": noise_model},
                                        "run_options":{"shots": shots}}
                    )
            else:
                estimator = Estimator(options={"run_options":{"shots": shots}})
            result = estimator.run([(qc, observable)]).result()

            C[i] -= 1 / 2 * H.coeffs[pauli_string].real * result[0].data.evs
    return np.array(C)


def pauli_measure(qc: QuantumCircuit, pauli_string: str) -> None:
    """Measure the given pauli string on the provided quantum circuit using a hadamard test.

    Args:
        qc (QuantumCircuit): The quantum circuit ansatz being constructed.
        pauli_string (str): The pauli string to be measured as a string.
    """
    N = len(pauli_string)
    for i in range(len(pauli_string)):  # Measure Pauli Strings
        if str(pauli_string[i]) == "X":
            qc.cx(N, i)
        if str(pauli_string[i]) == "Y":
            qc.cy(N, i)
        if str(pauli_string[i]) == "Z":
            qc.cz(N, i)


def Construct_Ansatz(
    init_circ: QuantumCircuit, params: npt.NDArray[np.float64], N: int
) -> QuantumCircuit:
    """Construct the full ansatz for use in measuring observables.

    Args:
        init_circ (QuantumCircuit): A qiskit circuit constructing the initial state of the system.
        params (numpy.array): A numpy vector containing the values of the parameters of the ansatz at a specific time.
        N (int): The number of qubits.

    Returns:
        QuantumCircuit: The full ansatz as a qiskit.QuantumCircuit.
    """
    qc = QuantumCircuit(N, 0)
    qc = qc.compose(init_circ, [k for k in range(N)])

    ansatz = QuantumCircuit(N, 0)
    for parameter in range(len(params)):  # Apply parameterized gates
        apply_param(params, parameter, ansatz, N)

    qc = qc.compose(ansatz, [k for k in range(N)])
    return qc


def ansatz_energy(
    init_circ: QuantumCircuit,
    params: npt.NDArray[np.float64],
    H: SparsePauliOp,
    shots: int = 2**14,
    noisy: bool = False,
) -> Tuple[float, float]:
    """Measure the energy of the ansatz.

    Args:
        init_circ (QuantumCircuit): A qiskit circuit constructing the initial state of the system.
        params (numpy.array): A numpy vector containing the values of the parameters of the ansatz at a specific time.
        H (SparsePauliOp): The Hamiltonian.
        shots (_type_, optional): The number of shots to be used to measure the energy. Defaults to 2**14.
        noisy (bool, optional): A boolean used to turn on and off the Fake-Sherbrooke qiskit noisy backend. Defaults to False.

    Returns:
        (float, float): Return (energy, variance) from the measured observables.
    """
    N = H.num_qubits

    if noisy:
        device_backend = FakeSherbrooke()
        coupling_map = device_backend.coupling_map
        noise_model = NoiseModel.from_backend(device_backend)
        basis_gates = noise_model.basis_gates
        estimator = Estimator(options={
                                "backend_options":{"noise_model": noise_model},
                                "run_options":{"shots": shots}}
                                )
    else:
        estimator = Estimator(options={"run_options":{"shots": shots}})
    qc = Construct_Ansatz(init_circ, params, N)
    result = estimator.run([(qc, H)]).result()
    return result[0].data.evs, result[0].data.stds


def VarQRTE(
    n_reps_ansatz: int,
    hamiltonian: SparsePauliOp,
    total_time: float = 1.0,
    timestep: float = 0.1,
    init_circ: Optional[QuantumCircuit] = None,
    shots: int = 2**10,
    noisy: bool = False,
) -> List[npt.NDArray[np.float64]]:
    """The Variational Quantum Real Time Evolution (VarQRTE) algorithm.  This uses quantum circuits to measure
        the elements of two objects, the A_ij matrix and the C_i vector.

    Args:
        n_reps_ansatz (int): The number of repetitions of the variational ansatz used to simulate Real-Time evolution.
        hamiltonian (SparsePauliOp): The Hamiltonian of the system.
        total_time (float, optional): A float to determine the total evolution time of the quantum system. Defaults to 1.0.
        timestep (float, optional): A float to determine the size of a single timestep. Defaults to 0.1.
        init_circ (QuantumCircuit, optional): A qiskit circuit constructing the initial state of the system.. Defaults to None.
        shots (int, optional): Number of shots to be used to measure observables. Defaults to 2**10.
        noisy (bool, optional): A boolean used to turn on and off the Fake-Sherbrooke qiskit noisy backend. Defaults to False.

    Returns:
        numpy.array: An array containing all the parameter values of the ansatz throughout its time evolution.
            These values can be put into Construct_Ansatz, or anstaz_energy to obtain observables of the system.
    """
    if init_circ is None:
        init_circ = QuantumCircuit(hamiltonian.num_qubits)

    initial_params = np.zeros(hamiltonian.num_qubits * (n_reps_ansatz + 1))
    num_timesteps = int(total_time / timestep)
    all_params = [np.copy(initial_params)]
    my_params = np.copy(initial_params)  # Reset Initial Parameters after each run
    for i in range(num_timesteps):
        print(f"Simulating Time={str(timestep*(i+1))}                      ", end="\r")
        theta_dot = np.array([0.0 for j in range(len(my_params))])
        A = Measure_A(
            init_circ, my_params, hamiltonian.num_qubits, shots=shots, noisy=noisy
        )
        C = Measure_C(
            init_circ,
            my_params,
            hamiltonian,
            hamiltonian.num_qubits,
            shots=shots,
            evolution_type="real",
            noisy=noisy,
        )

        # Approximately invert A using Truncated SVD
        u, s, v = np.linalg.svd(A)
        for j in range(len(s)):
            if s[j] < 1e-2:
                s[j] = 1e8
        t = np.diag(s**-1)
        A_inv = np.dot(v.transpose(), np.dot(t, u.transpose()))

        theta_dot = np.matmul(A_inv, C)

        my_params -= theta_dot * timestep
        all_params.append(np.copy(my_params))
    return all_params


def VarQITE(
    n_reps_ansatz: int,
    hamiltonian: SparsePauliOp,
    total_time: float,
    timestep: float,
    init_circ: Optional[QuantumCircuit] = None,
    shots: int = 2**10,
    noisy: bool = False,
) -> List[npt.NDArray[np.float64]]:
    """The Variational Quantum Imaginary Time Evolution (VarQITE) algorithm.  This uses quantum circuits to measure
        the elements of two objects, the A_ij matrix and the C_i vector.

    Args:
        n_reps_ansatz (int): The number of repetitions of the variational ansatz used to simulate Imaginary-Time evolution.
        hamiltonian (SparsePauliOp): The Hamiltonian of the system.
        total_time (float, optional): A float to determine the total evolution time of the quantum system. Defaults to 1.0.
        timestep (float, optional): A float to determine the size of a single timestep. Defaults to 0.1.
        init_circ (QuantumCircuit, optional): A qiskit circuit constructing the initial state of the system.. Defaults to None.
        shots (int, optional): Number of shots to be used to measure observables. Defaults to 2**10.
        noisy (bool, optional): A boolean used to turn on and off the Fake-Sherbrooke qiskit noisy backend. Defaults to False.

    Returns:
        numpy.array: An array containing all the parameter values of the ansatz throughout its time evolution.
            These values can be put into Construct_Ansatz, or anstaz_energy to obtain observables of the system.
    """
    if init_circ is None:
        init_circ = QuantumCircuit(hamiltonian.num_qubits)

    initial_params = np.zeros(hamiltonian.num_qubits * (n_reps_ansatz + 1))
    num_timesteps = int(total_time / timestep)
    all_params = [np.copy(initial_params)]

    my_params = np.copy(initial_params)  # Reset Initial Parameters after each run
    for i in range(num_timesteps):
        print(f"Timestep: {str(i*timestep)}                      ", end="\r")
        theta_dot = np.array([0.0 for j in range(len(my_params))])
        A = np.array(
            Measure_A(
                init_circ, my_params, hamiltonian.num_qubits, shots=shots, noisy=noisy
            )
        )
        C = np.array(
            Measure_C(
                init_circ,
                my_params,
                hamiltonian,
                hamiltonian.num_qubits,
                shots=shots,
                noisy=noisy,
                evolution_type="imaginary",
            )
        )

        # Approximately invert A using Truncated SVD
        u, s, v = np.linalg.svd(A)
        for j in range(len(s)):
            if s[j] < 1e-2:
                s[j] = 1e7
        t = np.diag(s**-1)
        A_inv = np.dot(v.transpose(), np.dot(t, u.transpose()))
        # A_inv=np.dot(v,np.dot(t,u.transpose()))

        theta_dot = np.matmul(A_inv, C)

        my_params += theta_dot * timestep
        all_params.append(np.copy(my_params))

        # print("Theta dot: "+str(np.sum(np.abs(theta_dot))))
        # print("(Energy ,Variance): "+str(ansatz_energy(init_circ, my_params[:], hamiltonian)))
        # print()
    return all_params
