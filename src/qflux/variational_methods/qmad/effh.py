import numpy as np
from scipy.linalg import expm, kron
 

def vectorize_comm(A):
    # Create an identity matrix with the same dimension as A
    iden = np.eye(A.shape[0])
    # Compute the vectorized commutator [A, .] as the Kronecker product 
    return kron(iden, A) - kron(A.T, iden)

class VectorizedEffectiveHamiltonian_class:
    def __init__(self, He, Ha):
        self.He = He
        self.Ha = Ha

def VectorizedEffectiveHamiltonian(H, gamma, lind):
    # Create an identity matrix with the same dimension as H
    iden = np.eye(H.shape[0])
    # Get the dimension of H
    d = H.shape[0]
    # Compute the vectorized commutator for the Hamiltonian H
    vec_H = vectorize_comm(H)
    # Initialize the result matrix with zeros (complex type)
    res = np.zeros((d**2, d**2), dtype=np.complex128)
    
    # Compute the conjugate of the Lindblad operator
    L_conj = lind.conj()
    L_dagger_L = L_conj.T @ lind
    # Compute the Lindblad contribution to the effective Hamiltonian
    res -= gamma * (kron(L_conj, lind) - (kron(iden, L_dagger_L) + kron(L_dagger_L.T, iden)) / 2)

    # Return an instance of the VectorizedEffectiveHamiltonian_class with vec_H and res
    return VectorizedEffectiveHamiltonian_class(vec_H, res)

class EffectiveHamiltonian_class:
    def __init__(self, He, Ha, Llist, LdL):
        self.He = He  # Hermitian part
        self.Ha = Ha  # Anti-Hermitian part
        self.Llist = Llist  # List of Lindblad operators
        self.LdL = LdL  # List of L†L

def EffectiveHamiltonian( mats, Llist):
    """
    Create an EffectiveHamiltonian object based on provided parameters.
    :param mats: List of matrices (Hamiltonian terms).
    :param Llist: List of lists of Lindblad operators.
    :return: An instance of EffectiveHamiltonian_class.
    """
    # Directly use the input matrices as the Hermitian part (assuming they are Hermitian)
    He = sum(mats)  # Sum of Hamiltonian terms as Hermitian part
    Ha = 0.0  # Initialize anti-Hermitian part
    LdL = []  # Initialize the list for Lindblad operator products

    for  LL in Llist:
        for L in LL:
            L_dagger_L = (L.conj().T @ L) 
            
            LdL.append(L_dagger_L)  # Append to LdL list
            Ha += L_dagger_L  # Sum for the anti-Hermitian part

    # Return the Effective Hamiltonian object
    return EffectiveHamiltonian_class(He, 0.5 * Ha, [L for LL in Llist for L in LL], LdL)

# Helper methods that align with the Julia logic
def herm(H, t):
    return H.He(t)  # Assuming He is callable with t

def antiherm(H):
    return H.Ha

def get_LdL(H):
    return H.LdL

def get_L(H):
    return H.Llist
