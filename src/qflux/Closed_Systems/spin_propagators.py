from qiskit.circuit.library import PauliEvolutionGate
# Trotter-Suzuki implementation for decomposition of exponentials
# of matrices
from qiskit.synthesis import SuzukiTrotter
from qiskit.quantum_info import SparsePauliOp
from qiskit import QuantumCircuit, QuantumRegister
import numpy as np
from itertools import groupby
import re



# -----------------------------------------------------------------
# Hamiltonian Functions
# -----------------------------------------------------------------
def get_hamiltonian_n_site_terms(n, coeff, n_qubits):
    '''
        Assembles each term in the Hamiltonian based on their Pauli string
        representation and multiplying by the respective coefficient.
    '''
    XX_coeff = coeff[0]
    YY_coeff = coeff[1]
    ZZ_coeff = coeff[2]
    Z_coeff = coeff[3]

    XX_term = SparsePauliOp(("I" * n + "XX" + "I" * (n_qubits - 2 - n)))
    XX_term *= XX_coeff
    YY_term = SparsePauliOp(("I" * n + "YY" + "I" * (n_qubits - 2 - n)))
    YY_term *= YY_coeff
    ZZ_term = SparsePauliOp(("I" * n + "ZZ" + "I" * (n_qubits - 2 - n)))
    ZZ_term *= ZZ_coeff
    Z_term = SparsePauliOp(("I" * n + "Z" + "I" * (n_qubits - 1 - n)))
    Z_term *= Z_coeff

    return (XX_term + YY_term + ZZ_term + Z_term)


def get_heisenberg_hamiltonian(n_qubits, coeff=None):
    r'''
        Takes an integer number corresponding to number of spins/qubits
        and a list of sublists containing the necessary coefficients
        to assemble the complete Hamiltonian:
        $$
            H = \sum _i ^N h_z Z_i
                + \sum _i ^{N-1} (h_xx X_i X_{i+1}
                                  + h_yy Y_i Y_{i+1}
                                  + h_zz Z_i Z_{i+1}
                                 )
        $$
        Each sublist contains the [XX, YY, ZZ, Z] coefficients in this order.
        The last sublist should have the same shape, but only the Z component
        is used.
        If no coefficient list is provided, all are set to 1.
    '''

    # Three qubits because for 2 we get H_O = 0
    assert n_qubits >= 3

    if coeff == None:
        'Setting default values for the coefficients'
        coeff = [[1.0, 1.0, 1.0, 1.0] for i in range(n_qubits)]

    # Even terms of the Hamiltonian
    # (summing over individual pair-wise elements)
    H_E = sum((get_hamiltonian_n_site_terms(i, coeff[i], n_qubits)
               for i in range(0, n_qubits-1, 2)))

    # Odd terms of the Hamiltonian
    # (summing over individual pair-wise elements)
    H_O = sum((get_hamiltonian_n_site_terms(i, coeff[i], n_qubits)
               for i in range(1, n_qubits-1, 2)))

    # adding final Z term at the Nth site
    final_term = SparsePauliOp("I" * (n_qubits - 1) + "Z")
    final_term *= coeff[n_qubits-1][3]
    if (n_qubits % 2) == 0:
        H_E += final_term
    else:
        H_O += final_term

    # Returns the list of the two sets of terms
    return [H_E, H_O]


def get_time_evolution_operator(num_qubits, tau, trotter_steps, coeff=None):
    '''
        Given a number of qubits, generates the corresponding time-evolution for
        the Ising model with the same number of sites.

        Input:
            num_qubits (int): number of qubits, which should be equal to the
                number of spins in the chain
            evo_time (float): time parameter in time-evolution operator
            trotter_steps (int): number of time steps for the Suzuki-Trotter
                decomposition
            coeff (list of lists): parameters for each term in the Hamiltonian
                for each site ie ([[XX0, YY0, ZZ0, Z0], [XX1, YY1, ZZ1, Z1], ...])
        Returns:
            evo_op.definition: Trotterized time-evolution operator
    '''
    # Constructing the Hamiltonian here;
    # heisenberg_hamiltonian = [H_E, H_O]
    heisenberg_hamiltonian = get_heisenberg_hamiltonian(num_qubits,
                                                        coeff)

    # e^ (-i*H*evo_time), with Trotter decomposition
    # exp[(i * evo_time)*(IIIIXXIIII + IIIIYYIIII + IIIIZZIIII + IIIIZIIIII)]
    evo_op = PauliEvolutionGate(heisenberg_hamiltonian, tau,
                                synthesis=SuzukiTrotter(order=2,
                                    reps=trotter_steps))
    # The Trotter order=2 applies one set of the operators for
    # half a timestep, then the other set for a full timestep,
    # then the first step for another half a step note that reps
    # includes the number of repetitions of the Trotterized
    # operator higher number means more repetitions, and thus
    # allowing larger timestep
    return evo_op.definition


def find_string_pattern(pattern, string):
    match_list = []
    for m in re.finditer(pattern, string):
        match_list.append(m.start())
    return match_list


# efficient propagator for pauli hamiltonians
def sort_Pauli_by_symmetry(ham):
    '''
        Separates a qiskit PauliOp object terms into 1 and 2-qubit
        operators. Furthermore, 2-qubit operators are separated according
        to the parity of the index first non-identity operation.
    '''
    one_qubit_terms = []
    two_qubit_terms = []
    # separating the one-qubit from two-qubit terms
    for term in ham:
        matches = find_string_pattern('X|Y|Z', str(term.paulis[0]))
        pauli_string = term.paulis[0]
        coeff = np.real(term.coeffs[0])
        str_tag = pauli_string.to_label().replace('I', '')
        if len(matches) == 2:
            two_qubit_terms.append((pauli_string, coeff, matches, str_tag))
        elif len(matches) == 1:
            one_qubit_terms.append((pauli_string, coeff, matches, str_tag))

    # sorting the two-qubit terms according to index on which they act
    two_qubit_terms = sorted(two_qubit_terms, key=lambda x: x[2])
    # separating the even from the odd two-qubit terms
    even_two_qubit_terms = list(filter(lambda x: not x[2][0]%2, two_qubit_terms))
    odd_two_qubit_terms = list(filter(lambda x: x[2][0]%2, two_qubit_terms))

    even_two_qubit_terms = [list(v) for i, v in groupby(even_two_qubit_terms, lambda x: x[2][0])]
    odd_two_qubit_terms = [list(v) for i, v in groupby(odd_two_qubit_terms, lambda x: x[2][0])]

    return one_qubit_terms, even_two_qubit_terms, odd_two_qubit_terms


def generate_circ_pattern_1qubit(circ, term, delta_t):
    '''
        General 1-qubit gate for exponential of product identity and
        a single pauli gate.
        Only a single rotation operation is required, with the angle
        being related to the exponential argument:

        R_P(coeff) = exp(-i * coeff * P / 2)

        Where P is the Pauli gate and coeff encompasses the constant
        coefficient term
    '''
    coeff = 2 * term[1] * delta_t
    if term[3] == 'X':
        circ.rx(coeff, term[2])
    elif term[3] == 'Y':
        circ.ry(coeff, term[2])
    elif term[3] == 'Z':
        circ.rz(coeff, term[2])

    return circ


def generate_circ_pattern_2qubit(circ, term, delta_t):
    r'''
        General 2-qubit gate for exponential of Paulis. This is the
        optimal decomposition, based on a component of a U(4) operator.
        (see )

        The circuit structure is as follows:

        - ---- I ---- C - Rz(o) - X --- I --- C - Rz(pi/2) -
        - Rz(-pi/2) - X - Ry(p) - C - Ry(l) - X ---- I -----

        Where CX represent CNOT operations, R are rotation gates with angles,
        and I is the identity matrix. The angles are parameterized as follows:

        $ o = \theta = (\pi/2 - A) $
        $ p = \phi = (A - \pi/2) $
        $ l = \lambda = (\pi/2 - A) $

        Where A is the exponential argument.
    '''
    # wires to which to apply the operation
    wires = term[0][2]

    # angles to parameterize the circuit,
    # based on exponential argument
    if any('XX' in sublist for sublist in term):
        g_phi = ( 2 * (-1) * term[0][1] * delta_t - np.pi / 2)
    else:
        g_phi = - np.pi / 2
    if any('YY' in sublist for sublist in term):
        g_lambda = (np.pi/2 - 2 * (-1) * term[1][1] * delta_t)
    else:
        g_lambda = np.pi/2
    if any('ZZ' in sublist for sublist in term):
        g_theta = (np.pi/2 - 2 * (-1) * term[2][1] * delta_t)
    else:
        g_theta = np.pi/2

    # circuit
    circ.rz(-np.pi/2, wires[1])
    circ.cx(wires[1], wires[0])
    circ.rz(g_theta, wires[0])
    circ.ry(g_phi, wires[1])
    circ.cx(wires[0], wires[1])
    circ.ry(g_lambda, wires[1])
    circ.cx(wires[1], wires[0])
    circ.rz(np.pi/2, wires[0])
    return circ


def get_manual_Trotter(num_q, pauli_ops, timestep, n_trotter=1,
                       trotter_type='basic', reverse_bits=True):
    # sorts the Pauli strings according to qubit number they affect and symmetry
    one_q, even_two_q, odd_two_q = sort_Pauli_by_symmetry(pauli_ops)
    # scales the timestep according to the number of trotter steps
    timestep_even_two_q = timestep / n_trotter
    timestep_odd_two_q = timestep / n_trotter
    timestep_one_q = timestep / n_trotter
    # symmetric places 1/2 of one_q and odd_two_q before and after even_two_q
    if trotter_type == 'symmetric':
        timestep_odd_two_q /= 2
        timestep_one_q /= 2
    # constructs circuits for each segment of the operators
    qc_odd_two_q, qc_even_two_q, qc_one_q = QuantumCircuit(num_q), QuantumCircuit(num_q), QuantumCircuit(num_q)
    for i in even_two_q:
        qc_even_two_q = generate_circ_pattern_2qubit(qc_even_two_q, i, timestep_even_two_q)
    for i in odd_two_q:
        qc_odd_two_q = generate_circ_pattern_2qubit(qc_odd_two_q, i, timestep_odd_two_q)
    for i in one_q:
        qc_one_q = generate_circ_pattern_1qubit(qc_one_q, i, timestep_one_q)
    # assembles the circuit for Trotter decomposition of exponential
    qr = QuantumRegister(num_q)
    qc = QuantumCircuit(qr)
    if trotter_type == 'basic':
        qc = qc.compose(qc_even_two_q)
        qc = qc.compose(qc_odd_two_q)
        qc = qc.compose(qc_one_q)
    elif trotter_type == 'symmetric':
        qc = qc.compose(qc_one_q)
        qc = qc.compose(qc_odd_two_q)
        qc = qc.compose(qc_even_two_q)
        qc = qc.compose(qc_odd_two_q)
        qc = qc.compose(qc_one_q)
    # repeats the single_trotter circuit several times to match n_trotter
    for i in range(n_trotter-1):
        qc = qc.compose(qc)
    if reverse_bits:
        return qc.reverse_bits()
    else:
        return qc


if __name__ == '__main__':
    num_shots = 100
    num_q = 3
    evolution_timestep = 0.1
    n_trotter_steps = 1
    # XX YY ZZ, Z
    ham_coeffs = ([[0.75/2, 0.75/2, 0.0, 0.65]]
                    + [[0.5, 0.5, 0.0, 1.0]
                    for i in range(num_q-1)])
    time_evo_op = get_time_evolution_operator(
        num_qubits=num_q, tau=evolution_timestep,
        trotter_steps=n_trotter_steps, coeff=ham_coeffs)
    print(time_evo_op)

    spin_chain_hamiltonian = get_heisenberg_hamiltonian(num_q,
                                                        ham_coeffs)

    spin_chain_hamiltonian = sum(spin_chain_hamiltonian)
    print(get_manual_Trotter(num_q, spin_chain_hamiltonian,
                             0.1).draw())
    print(get_manual_Trotter(num_q, spin_chain_hamiltonian, 0.1,
                             n_trotter=2).draw())
    print(get_manual_Trotter(num_q, spin_chain_hamiltonian, 0.1,
                             trotter_type='symmetric').draw())
    print(get_manual_Trotter(num_q, spin_chain_hamiltonian, 0.1,
                             n_trotter=2,
                             trotter_type='symmetric').draw())
