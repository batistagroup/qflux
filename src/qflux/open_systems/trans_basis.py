import numpy as np
import scipy.fft as sfft

#the module in the current package
from . import params as pa
    
def x2k_wave(dx,psi):
    """
    transform the wavefunction from x space to k space
    dx: the interval of the x-space grid point
    psi: the wavefunction in the x-space
    """
    pre_fac = dx/(2*np.pi)**0.5
    psik = sfft.fft(psi)*pre_fac
    return psik

def k2x_wave(dx,psik):
    """
    transform the wavefunction from x space to k space
    dx: the interval of the x-space grid point
    psi: the wavefunction in the k-space
    """
    pre_fac = (2*np.pi)**0.5/dx
    psix  = sfft.ifft(psik.copy())*pre_fac
    return psix

def nested_kronecker_product(pauli_str):
    '''
    Handles Kronecker Products for list (i.e., pauli_str = 'ZZZ' will evaluate Z Z Z).
    Given string 'pauli_str' this evaluates the kronecker product of all elements.
    '''

    # Define a dictionary with the four Pauli matrices:
    pms = {'I': pa.I,'X': pa.X,'Y': pa.Y,'Z': pa.Z}

    result = 1
    for i in range(len(pauli_str)):
     result = np.kron(result,pms[pauli_str[i]])
    return result


def pauli_to_ham(pauli_dict, Nqb):
    ''' 
        Function that Assembles the Hamiltonian based on their Pauli string
        
        pauli_dict: A dictionary that contain all the pauli string and the value of the Hamiltonian
            (the key is pauli string and the value is coefficient)
        
        Nqb: number of qubits, should match the length of the pauli string
        
        return: Hamiltonian matrix 
    '''
        
    Hmat = np.zeros((2**Nqb,2**Nqb),dtype=np.complex128)
    for key in pauli_dict:
        Hmat += pauli_dict[key]*nested_kronecker_product(key)

    return Hmat

def ham_to_pauli(Ham_arr, Nqb, tol=1E-5):
    '''
    Function that decomposes `Ham_arr` into a sum of Pauli strings.
    result: a dictionary with the key is pauli string and the value is coefficient
    '''
    import itertools

    pauli_keys = ['I','X','Y','Z'] # Keys of the dictionary

    if(2**Nqb != Ham_arr.shape[0]):
        print('Nqb and Matrix size not matched!')

    # Make all possible tensor products of Pauli matrices sigma
    sigma_combinations = list(itertools.product(pauli_keys, repeat=Nqb))

    result = {} # Initialize an empty dictionary to the results
    for ii in range(len(sigma_combinations)):
        pauli_str = ''.join(sigma_combinations[ii])

        # Evaluate the Kronecker product of the matrix array
        tmp_p_matrix = nested_kronecker_product(pauli_str)

        # Compute the coefficient for each Pauli string
        a_coeff = (1/(2**Nqb)) * np.trace(tmp_p_matrix @ Ham_arr)

        # If the coefficient is non-zero, we want to use it!
        if abs(a_coeff) > tol:
            result[pauli_str] = a_coeff.real

    return result

def trans_basis(operator,nbasis,psi_newbasis):
    """
    operator: the operator matrix in the old basis
    nbasis: the truncation in the new basis
    psi_newbasis: the new basis wave function expressed in the old basis 
               psi_newbasis[i,j] = <i|psi_newbasis|j>, with i is old basis, j new basis
    """
    
    operator_new = np.zeros((nbasis,nbasis),dtype=np.complex128)
    
    for i in range(nbasis):
      for j in range(nbasis):
        operator_new[i,j] = np.dot(np.dot(psi_newbasis[:,i].conj(),operator),psi_newbasis[:,j])

    return operator_new

def trans_basis_diag(diag_opr,nbasis,psi_newbasis):
    """
    diag_opr: the array represent the diagnoal operator in the old basis
    nbasis: the truncation in the new basis
    psi_newbasis: the new basis wave function expressed in the old basis 
               psi_newbasis[i,j] = <i|psi_newbasis|j>, with i is old basis, j new basis
    """
    
    operator_new = np.zeros((nbasis,nbasis),dtype=np.complex128)
    
    for i in range(nbasis):
      for j in range(nbasis):
        operator_new[i,j] = np.dot(np.multiply(psi_newbasis[:,i].conj(),diag_opr),psi_newbasis[:,j])

    return operator_new