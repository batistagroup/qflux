# Utilities
# Plotting utilities, conversion factors, etc.
import numpy as np
import qiskit_aer
from qiskit.compiler import transpile
from qiskit_ibm_runtime import Sampler
from qiskit.quantum_info import SparsePauliOp
import itertools

# Conversion Factors:
def convert_au_to_eV(input_val):
    au2ev = 27.21138602
    return(au2ev * input_val)


def convert_eV_to_au(input_val):
    au2ev = 27.21138602
    ev2au = 1/au2ev
    return(ev2au * input_val)


def convert_bohr_to_au(input_val):
    bohr2au = 0.52917721092
    return(bohr2au * input_val)


def convert_au_to_bohr(input_val):
    bohr2au = 0.52917721092
    au2bohr = 1/bohr2au
    return(input_val * au2bohr)


def convert_fs_to_au(input_val):
    fs2au = 41.3414
    return(fs2au * input_val)


def convert_au_to_fs(input_val):
    fs2au = 41.3414
    au2fs = 1/fs2au
    return(au2fs * input_val)


def get_proton_mass():
    proton_mass = 1836.15267343 # proton-electron mass ratio
    return(proton_mass)


def execute(QCircuit, backend=None, shots=None, real_backend=False):
    '''
        Function to replace the now-deprecated Qiskit
        `QuantumCircuit.execute()` method.

        Input:
          - `QCircuit`: qiskit.QuantumCircuit object
          - `Backend`: qiskit.Backend instance
          - `shots`: int specifying the number of shots
          - `real_backend`: bool specifying whether the provided backend is
            a real device (True) or not (False)
    '''
    if shots:
        n_shots = shots
    else:
        n_shots = 1024 # Use the qiskit default if not specified

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


# Calculation of Expectation Value:
def calculate_expectation_values(dynamics_results, observable_grid, do_FFT=False, dx=None):
    '''
    Function to calculate the time-dependent expectation value of an observable O defined on a grid.
    Inputs:

        - `dynamics_results`: np.ndarray of wavefunctions/propagated states with shape: (n_steps, nx)
        - `observable_grid`: np.array of observable
    '''
    if dx:
        d_observable = dx
    else:
        d_observable = observable_grid[1] - observable_grid[0]
    if do_FFT:
        psi_list = np.fft.fft(dynamics_results, axis=1, norm='ortho')
    else:
        psi_list = dynamics_results
    # Compute the expectation value.
    expectation  = np.real(np.sum(psi_list.conj()*observable_grid*psi_list*d_observable, axis=1))

    return(expectation)

# Pauli Decomposition Utilities

def vec_query(arr, my_dict):
    '''
    This function vectorizes dictionary querying, allowing us to query `my_dict` with a np.array `arr` of keys.
    '''

    return np.vectorize(my_dict.__getitem__, otypes=[tuple])(arr)


def nested_kronecker_product(a):
    '''
    Handles Kronecker Products for list (i.e.,  a = [Z, Z, Z] will evaluate Z ⊗ Z ⊗ Z)
    '''
    if len(a) == 2:
            return np.kron(a[0],a[1])
    else:
        return np.kron(a[0], nested_kronecker_product(a[1:]))


def Hilbert_Schmidt(mat1, mat2):
    'Return the Hilbert-Schmidt Inner Product of two matrices.'
    return np.trace(mat1.conj().T @ mat2)


def decompose(Ham_arr, tol=12, subset = None):
    '''
    Function that takes an input matrix `H` and decomposes into a sum of tensor products of Pauli Matrices.
        - The expectation is that `H` will have a shape of 2^n, where n is the number of qubits needed to represent
        the matrix in this form. n determines the number of terms in the tensor product composing the Pauli strings.
        - Ex: If H has shape of (16, 16), n = log2(H) = 4 qubits.
         The Pauli strings will then take the form of ['IIII', 'ZIII', 'ZZII', etc.]
    Has the option to toggle verbose output, which prints the terms of the pauli sum.
    '''

    X = np.asarray([[0,1],[1,0]])
    Y = np.asarray([[0,complex(0,-1)],[complex(0,1),0]])
    Z = np.asarray([[1,0],[0,-1]])
    I = Z@Z
    # Define a dictionary with the four Pauli matrices:
    global_pms = {'I': I,'X': X,'Y': Y,'Z': Z}

    # The next few lines allow for the function use a reduced basis.
    pm_keys = []
    pm_vals = []
    if subset:
        for tmp_key in subset:
            pm_keys.append(tmp_key)
            pm_vals.append(global_pms[tmp_key])
        pms = dict(zip(pm_keys, pm_vals))
    else:
        pms = global_pms

    pauli_keys = list(pms.keys()) # Keys of the dictionary

    nqb = int(np.log2(Ham_arr.shape[0])) # Determine the # of qubits needed

    # Make all possible tensor products of Pauli matrices sigma
    sigma_combinations = list(itertools.product(pauli_keys, repeat=nqb))

    output_string = '' # Initialize an empty string to which we can add our terms
    for ii in range(len(sigma_combinations)):
        pauli_str = ''.join(sigma_combinations[ii])

        # Convert the Pauli string into a list of matrices
        tmp_mat_list = vec_query(np.array(list(sigma_combinations[ii])), pms)

        # Evaluate the Kronecker product of the matrix array
        tmp_p_matrix = nested_kronecker_product(tmp_mat_list)

        # Compute the coefficient for each Pauli string
        a_coeff = (1/(2**nqb)) * Hilbert_Schmidt(tmp_p_matrix, Ham_arr)

        # If the coefficient is non-zero, we want to use it!
        if abs(a_coeff) > 10**(-tol):
            output_string += str(np.round(a_coeff.real, tol))+'*'+pauli_str
            output_string += '+' # Add a plus sign for the next term!

    return output_string[:-1] # To ignore that extra plus sign


def build_pauli_dict(decomposed_operator):
    '''
    Build a Pauli dictionary from the output of the `decompose` function above.
    This is done by converting the unique Pauli Strings ['IIII', 'IIIZ', etc.] to keys,
    where the values of the dictionary are the numerical coefficients.
    '''
    single_terms = decomposed_operator.split('+')
    pauli_dict_out = {}
    for term in single_terms:
        coeff, pauli_str = term.split('*')
        pauli_dict_out[pauli_str] = float(coeff)
    return pauli_dict_out


def pauli_strings_2_pauli_sum(operator):
    '''
    Function to convert the output of the `decompose` function above into a Qiskit-recognized PauliSumOp.
    The string representing the Hamiltonian as a sum of product of Paulis is converted into a dictionary where
    the keys are the Pauli Strings (ex: 'IIII' or 'IXYZ') and the values are the coefficients.
    '''
    tmp_pauli_dict = build_pauli_dict(operator)
    # Convert the dict to a list of tuples of the form [('Pauli String', float_coeff), ...]
    tmp_pauli_sum = SparsePauliOp(data=list(tmp_pauli_dict.keys()), coeffs=list(tmp_pauli_dict.values()))
    return tmp_pauli_sum
