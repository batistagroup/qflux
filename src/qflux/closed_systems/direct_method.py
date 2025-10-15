import numpy as np
from qiskit import QuantumCircuit, QuantumRegister

def exp_all_z(circuit, quantum_register, pauli_indexes, control_qubit=None, t=1):
    r"""
    Implements \( e^{-i t Z \otimes \cdots \otimes Z} \) on specified qubits.
    
    Args:
        circuit (QuantumCircuit): The circuit to modify.
        quantum_register (QuantumRegister): Register containing target qubits.
        pauli_indexes (list): Indices of qubits where \( Z \) acts.
        control_qubit (Qubit, optional): Optional control qubit for conditional application.
        t (float): Evolution time.
    
    Returns:
        QuantumCircuit: Updated circuit with the operation applied.
    """
    if control_qubit and control_qubit.register not in circuit.qregs:
        circuit.add_register(control_qubit.register)

    if not pauli_indexes:
        if control_qubit:
            circuit.p(t, control_qubit)  # Phase gate
        return circuit

    # Parity computation
    for i in range(len(pauli_indexes) - 1):
        circuit.cx(quantum_register[pauli_indexes[i]], quantum_register[pauli_indexes[i + 1]])

    # Apply phase rotation
    target = quantum_register[pauli_indexes[-1]]
    angle = -2 * t
    if control_qubit:
        circuit.crz(angle, control_qubit, target)
    else:
        circuit.rz(angle, target)

    # Uncompute parity
    for i in reversed(range(len(pauli_indexes) - 1)):
        circuit.cx(quantum_register[pauli_indexes[i]], quantum_register[pauli_indexes[i + 1]])

    return circuit


def exp_pauli(pauli, quantum_register, control_qubit=None, t=1):
    r"""
    Implements \( e^{-i t P} \) for a Pauli string \( P \).

    Args:
        pauli (str): Pauli string (e.g., "XIZY").
        quantum_register (QuantumRegister): Target register.
        control_qubit (Qubit, optional): Optional control qubit.
        t (float): Evolution time.

    Returns:
        QuantumCircuit: Circuit implementing the Pauli evolution.
    """
    if len(pauli) != len(quantum_register):
        raise ValueError("Pauli string length must match register size.")

    pauli_indexes = []
    pre_circuit = QuantumCircuit(quantum_register)

    for i, op in enumerate(pauli):
        if op == 'I':
            continue
        elif op == 'X':
            pre_circuit.h(i)
            pauli_indexes.append(i)
        elif op == 'Y':
            pre_circuit.rx(np.pi/2, i)
            pauli_indexes.append(i)
        elif op == 'Z':
            pauli_indexes.append(i)
        else:
            raise ValueError(f"Invalid Pauli operator '{op}' at position {i}.")

    circuit = QuantumCircuit(quantum_register)
    circuit.compose(pre_circuit, inplace=True)
    circuit = exp_all_z(circuit, quantum_register, pauli_indexes, control_qubit, t)
    circuit.compose(pre_circuit.inverse(), inplace=True)
    return circuit


def hamiltonian_simulation(hamiltonian, quantum_register=None, control_qubit=None, t=1, trotter_number=1):
    r"""
    Implements \( e^{-i H t} \) using first-order Trotterization.
    
    Args:
        hamiltonian (dict): Pauli terms with coefficients (e.g., {"ZZ": 0.5, "XX": 0.3}).
        quantum_register (QuantumRegister, optional): Target register.
        control_qubit (Qubit, optional): Optional control qubit.
        t (float): Simulation time.
        trotter_number (int): Number of Trotter steps.
    
    Returns:
        QuantumCircuit: Trotterized Hamiltonian evolution circuit.
    """
    if not hamiltonian:
        raise ValueError("Hamiltonian must contain at least one term.")
    
    n_qubits = len(next(iter(hamiltonian)))
    if quantum_register is None:
        quantum_register = QuantumRegister(n_qubits)
    
    delta_t = t / trotter_number
    circuit = QuantumCircuit(quantum_register)
    
    for pauli_str, coeff in hamiltonian.items():
        term_circuit = exp_pauli(pauli_str, quantum_register, control_qubit, coeff * delta_t)
        circuit.compose(term_circuit, inplace=True)
    
    full_circuit = QuantumCircuit(quantum_register)
    for _ in range(trotter_number):
        full_circuit.compose(circuit, inplace=True)
    
    return full_circuit


