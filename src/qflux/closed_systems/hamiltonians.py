from qiskit.quantum_info import SparsePauliOp


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
                + \sum _i ^{N-1} (h_xx X_iX_{i+1}
                                  + h_yy Y_iY_{i+1}
                                  + h_zz Z_iZ_{i+1}
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


if __name__ == '__main__':
    num_q = 3
    # XX YY ZZ, Z
    ham_coeffs = ([[0.75/2, 0.75/2, 0.0, 0.65]]
                  + [[0.5, 0.5, 0.0, 1.0]
                  for i in range(num_q-1)])

    spin_chain_hamiltonian = get_heisenberg_hamiltonian(num_q,
                                                        ham_coeffs)
    print('Hamiltonian separated into even and odd components:')
    print(spin_chain_hamiltonian)
    print('Hamiltonian combining even and odd components:')
    print(sum(spin_chain_hamiltonian))
