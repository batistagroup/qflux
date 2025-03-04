import numpy as np
import scipy.linalg as LA

from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit
from qiskit.quantum_info.operators import Operator

from . import walsh_gray_optimization as wo


def scale_array(array,scale=1.1):
    """
    renormalize an array to make sure it is a contraction
    """
    # Normalization factor, divide by martix's norm to ensure contraction
    # may divide a larger scaling factor controlled by the number (scale)
    norm = LA.norm(array,2)*scale
    array_new = array/norm
    
    return array_new, norm

def dilate_Sz_Nagy(array):
    """
    dilate the non-unitary array (should be a contraction) to a unitary matrix
    array: ndarray of N*N
    """
    ident = np.eye(array.shape[0])
      
    # Calculate the conjugate transpose of the G propagator
    fcon = (array.conjugate()).T
      
    # Calculate the defect matrix for dilation
    fdef = LA.sqrtm(ident - np.dot(fcon, array))
      
    # Calculate the defect matrix for the conjugate of the G propagator
    fcondef = LA.sqrtm(ident - np.dot(array, fcon))
      
    # Dilate the G propagator to create a unitary operator
    array_dilated = np.block([[array, fcondef], [fdef, -fcon]])
    
    return array_dilated
          
#generate the quantum gate matrices for SVD-dilation
def dilate_SVD(array):
    """
    dilate the non-unitary array to unitary matrices using SVD technique by Schlimgen et al
    (Phys. Rev. A 2022, 106, 022414.)
    
    array: ndarray of N*N (should be a contraction)
    """
    
    #get the dimension of the array
    Nvec = array.shape[0]
    
    #Performing SVD to the array
    U1,S1,V1 = LA.svd(array)

    Mzero = np.zeros((Nvec,Nvec),dtype=np.complex128)
    fk=np.zeros(2*Nvec,dtype=np.float64)
    
    Sig_p = np.zeros(Nvec,dtype=np.complex128)
    Sig_m = np.zeros(Nvec,dtype=np.complex128)
    for i in range(len(S1)):
        
        Sig_p[i] = S1[i]+1j*np.sqrt((1-S1[i]**2))
        Sig_m[i] = S1[i]-1j*np.sqrt((1-S1[i]**2))
    
        #here U_Sigma = e^{i f_k}
        fk[i] = (-1j*np.log(Sig_p[i])).real
        fk[Nvec+i] = (-1j*np.log(Sig_m[i])).real
    
    SG0 = np.block([[np.diag(Sig_p),Mzero],\
          [Mzero,np.diag(Sig_m)]])

    return U1, V1, SG0, fk


def cons_SVD_cirq(Nqb, array, ini_vec, Iswalsh = True):
    """
    Construct the SVD-dilation circuit
    Nqb: number of qubits
    array: the non-unitary array (should be a contraction)
    ini_vec: initial qubit state vector (in the dilation space)
    Iswalsh: If True, then combine with Walsh operator representation to reduce the circuit depth
    """
    
    qc = QuantumCircuit(Nqb,Nqb)
    qc.initialize(ini_vec,range(0,Nqb))
    
    U_matu,U_matv,S_mat0,fk = dilate_SVD(array)
    
    qc.append(Operator(U_matv),range(0,Nqb-1))
    qc.h(Nqb-1)

    if(Iswalsh):
        #the walsh coeff
        arr_a = wo.walsh_coef(fk,Nqb)
        
        Ulist_diag0 = wo.cirq_list_walsh(arr_a,Nqb,1E-5)
        Ulist_diag = wo.optimize(Ulist_diag0)
        qc_diag = wo.cirq_from_U(Ulist_diag,Nqb)
        
        qc.append(qc_diag.to_gate(),range(Nqb))

    else:
        qc.append(Operator(S_mat0),range(Nqb))
      
    qc.append(Operator(U_matu),range(0,Nqb-1))
    qc.h(Nqb-1)

    return qc

def cons_SzNagy_cirq(Nqb, array, ini_vec):
    """
    Construct the Sz-Nagy-dilation circuit
    Nqb: number of qubits
    array: the non-unitary array (should be a contraction)
    ini_vec: initial qubit state vector (in the dilation space)
    """
    
    qr = QuantumRegister(Nqb)  # Create a quantum register
    cr = ClassicalRegister(Nqb)  # Create a classical register to store measurement results
    qc = QuantumCircuit(qr, cr)  # Combine the quantum and classical registers to create the quantum circuit
    
    # Initialize the quantum circuit with the initial state
    qc.initialize(ini_vec, qr)
    
    # Create a custom unitary operator with the dilated propagator
    U_dil = dilate_Sz_Nagy(array)
    U_dil_op = Operator(U_dil)
    
    # Apply the unitary operator to the quantum circuit's qubits
    qc.unitary(U_dil_op, qr)

    return qc

def construct_circuit(Nqb, array,statevec,method = 'Sz-Nagy',Isscale = True):
    """
    construct the quantum circuit
    
    method: specify the dilation method, can be 'Sz-Nagy' or 'SVD' or 'SVD-Walsh'
        'Sz-Nagy': using Sz-Nagy method do dilation
        'SVD': using SVD-dilation method 
        'SVD-Walsh': SVD-dilation combine with Walsh operator representation to reduce the circuit depth
    
    Nqb: number of qubits in the original space, should match the dimension of array (N=2^Nqb)
    array: ndarray of N*N
    statevec: ndarray of N*1 (initial statevector)
    
    Isscale: 
        if Ture, renormalize the array to make sure it is a contraction
        if False, do not renormalize. 
    """   
    
    #the number of qubits for the original space
    if(2**Nqb != array.shape[0]):
        print('error, array dimension not matched!')
    
    #state vector in the dilated space
    statevec_dil = np.concatenate((statevec, np.zeros_like(statevec)))
    
    #first scale the array to make it a contraction
    if(Isscale):
        array_new, normfac = scale_array(array)
    else:
        array_new = array
    
    if(method == 'Sz-Nagy'):
        qc = cons_SzNagy_cirq(Nqb+1, array_new, statevec_dil)
    elif(method == 'SVD'):
        Iswalsh = False
        qc = cons_SVD_cirq(Nqb+1, array_new, statevec_dil, Iswalsh)
    elif(method == 'SVD-Walsh'):
        Iswalsh = True
        qc = cons_SVD_cirq(Nqb+1, array_new, statevec_dil, Iswalsh)
    else:
        print('method error in construct circuit')
    
    if(Isscale):
        return qc,normfac
    else:
        return qc


