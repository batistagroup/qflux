"""
This section uses some MPSQD functions to perform TT-TFD calculations for the Spin-Boson model dynamics. 
Created by Xiaohan Dan
"""

import numpy as np
from . import params as pa
import sys

from mpsqd.utils import MPS,add_tensor, MPS2MPO, calc_overlap
from .tdvp import tdvp1site
import time

 
def initial(istate):
        
    """
    initialize the state in the tensor train format For TT-TFD calculation 
    
    """
      
    #initial state 
    #Build initial ground state at spin-Up state
    su = np.zeros((1,pa.DOF_E,pa.MAX_TT_RANK),dtype=np.complex128)
    sd = np.zeros((1,pa.DOF_E,pa.MAX_TT_RANK),dtype=np.complex128)
    
    su[0,:,0] = pa.spin_up
    sd[0,:,0] = pa.spin_down
    
    e1 = np.sqrt(0.5) * (su + sd)
    e2 = np.sqrt(0.5) * (su + 1j * sd)
    
    # initial MPS wavepacket 
    nbarr = np.ones((1+2*pa.DOF_N))*pa.occ
    nbarr[0] = pa.DOF_E
    y0 = MPS(1+2*pa.DOF_N,nb=nbarr)  
    
    if(istate==0):
        y0.nodes.append(su)
    elif(istate==1):
        y0.nodes.append(e1)
    elif(istate==2):
        y0.nodes.append(e2)
    elif(istate==3):
        y0.nodes.append(sd)
       
    gs = np.zeros((pa.MAX_TT_RANK,pa.occ,pa.MAX_TT_RANK),dtype=np.complex128)
    gs[0,0,0] = 1.
    for k in range(2 * pa.DOF_N-1): # double space formation
        y0.nodes.append(gs)
    gs = np.zeros((pa.MAX_TT_RANK,pa.occ,1),dtype=np.complex128)
    gs[0,0,0] = 1.
    y0.nodes.append(gs)

    return y0


def construct_Hamil(eps=1E-14):
    """
    construct -iH for TT-TFD calculation, H is the effective Hamiltonian 
    Spin-Boson model, Ohmic spectral density
    
    eps: the truncation precision for effective Hamiltonian
    
    Return: an MPO object. 
    """
    
    om = pa.OMEGA_C / pa.DOF_N * (1 - np.exp(-pa.OMEGA_MAX/pa.OMEGA_C))

    # initialize arrays for parameters
    freq = np.zeros((pa.DOF_N)) # frequency
    ck = np.zeros((pa.DOF_N))   # linear electron-phonon coupling constant
    gk = np.zeros((pa.DOF_N))   # ck in occupation number representation
    thetak = np.zeros((pa.DOF_N)) # temperature-dependent mixing parameter in TFD
    sinhthetak = np.zeros((pa.DOF_N)) # sinh(theta)
    coshthetak = np.zeros((pa.DOF_N)) # cosh(theta)
    for i in range(pa.DOF_N):
        freq[i] = -pa.OMEGA_C * np.log(1-(i+1) * om/(pa.OMEGA_C)) # Ohmic frequency
        ck[i] = np.sqrt(pa.XI * om) * freq[i] #Ohmic coupling constant
        gk[i] = -ck[i] / np.sqrt(2 * freq[i]) #Transfer ck to occ. num. representation
    
        thetak[i] = np.arctanh(np.exp(-pa.BETA * freq[i]/2)) #theta, defined for harmonic models
        sinhthetak[i] = np.sinh(thetak[i]) #sinh(theta)
        coshthetak[i] = np.cosh(thetak[i]) #cosh(theta)
  
    #MPO of identity for all nuclear DOFs
    ident_nuclear = tt_eye(pa.DOF_N*2, pa.occ)
  
    # constructing Pauli operators
    px = np.array([[0.0,1],[1,0]],dtype=np.complex128)
    pz = np.array([[1.0,0],[0,-1]],dtype=np.complex128)
    # Build electronic site energy matrix
    He = pa.EPSILON  * pz + pa.GAMMA_DA  * px
    # TT-ize that energy matrix
    tt_He = tt_matrix(He)
    tt_He = tt_kron(tt_He, ident_nuclear)
  

  
    # Build number operator, corresponds to harmonic oscillator Hamiltonian
    numoc = np.diag(np.arange(0, pa.occ, 1,dtype=np.complex128))
    # Initiate the TT-ized number operator as a zero TT array with shape of occ^N
    tt_numoc = tt_eye(pa.DOF_N,pa.occ)
    for i in range(pa.DOF_N): tt_numoc.nodes[i] *= 0.0
    
    # Construct number operator as TT
    for k in range(pa.DOF_N):
        tmp0 = tt_matrix(numoc)
        tmp0.nodes[0] *= freq[k]
        if k == 0:
            tmp = tt_kron(tmp0, tt_eye(pa.DOF_N-1, pa.occ))
        elif 0 < k < pa.DOF_N-1:
            tmp = tt_kron(tt_eye(k-1,pa.occ), tmp0)
            tmp = tt_kron(tmp,tt_eye(pa.DOF_N - k,pa.occ))
        else:
            tmp = tt_kron(tt_eye(k, pa.occ),tmp0)
        tt_numoc = add_tensor(tt_numoc, tmp,small=eps)


    # Ensure correct dimensionality
    tt_Ie = tt_matrix(np.eye(2,dtype=np.complex128))
    tt_systemnumoc = tt_kron(tt_Ie, tt_numoc)
    tt_systemnumoc = tt_kron(tt_systemnumoc, tt_eye(pa.DOF_N,pa.occ))
    
    # create a duplicate of number operator for the ficticious system
    tt_tildenumoc = tt_kron(tt_Ie, tt_eye(pa.DOF_N,pa.occ))
    tt_tildenumoc = tt_kron(tt_tildenumoc, tt_numoc)
    
    thetak = np.zeros((pa.DOF_N)) #temperature-dependent mixing parameter in TFD
    sinhthetak = np.zeros((pa.DOF_N)) #sinh(theta)
    coshthetak = np.zeros((pa.DOF_N)) #cosh(theta)
    for i in range(pa.DOF_N):
        thetak[i] = np.arctanh(np.exp(-pa.BETA * freq[i]/2)) #theta, defined for harmonic models
        sinhthetak[i] = np.sinh(thetak[i]) #sinh(theta)
        coshthetak[i] = np.cosh(thetak[i]) #cosh(theta)
    
    #Build displacement operator, corresponds to x operator in real space
    eneroc = np.zeros((pa.occ, pa.occ),dtype=np.complex128)
    for i in range(pa.occ - 1):
        eneroc[i,i+1] = np.sqrt(i+1)
        eneroc[i+1,i] = eneroc[i,i+1]
    
    # initialize displacement operator
    tt_energy = tt_eye(pa.DOF_N,pa.occ)
    for i in range(pa.DOF_N): tt_energy.nodes[i] *= 0.0
    
    for k in range(pa.DOF_N):
        tmp0 = tt_matrix(eneroc)
        tmp0.nodes[0] *= gk[k] * coshthetak[k]
        if k == 0:
            # coshtheta takes account for energy flow from real to ficticious system
            # thus takes account for temperature effect
            tmp = tt_kron(tmp0, tt_eye(pa.DOF_N-1, pa.occ))
        elif 0 < k < pa.DOF_N - 1:
            tmp = tt_kron(tt_eye(k-1,pa.occ), tmp0)
            tmp = tt_kron(tmp,tt_eye(pa.DOF_N - k,pa.occ))
        else:
            tmp = tt_kron(tt_eye(k, pa.occ), tmp0)
        
        tt_energy = add_tensor(tt_energy, tmp, small=eps)

    tt_systemenergy = tt_kron(tt_matrix(pz), tt_energy)
    tt_systemenergy = tt_kron(tt_systemenergy, tt_eye(pa.DOF_N, pa.occ))
  
  
    # initialize displacement operator
    tt_tilenergy = tt_eye(pa.DOF_N,pa.occ)
    for i in range(pa.DOF_N): tt_tilenergy.nodes[i] *= 0.0
    
    for k in range(pa.DOF_N):
        tmp0 = tt_matrix(eneroc)
        tmp0.nodes[0] *= gk[k] * sinhthetak[k]
        if k == 0:
            tmp = tt_kron(tmp0 , tt_eye(pa.DOF_N-1, pa.occ))
        elif 0 < k < pa.DOF_N - 1:
            tmp = tt_kron(tt_eye(k-1,pa.occ), tmp0)
            tmp = tt_kron(tmp, tt_eye(pa.DOF_N - k,pa.occ))
        else:
            tmp = tt_kron(tt_eye(k, pa.occ), tmp0)
        tt_tilenergy = add_tensor(tt_tilenergy, tmp,small=eps) 

    tt_tildeenergy = tt_kron(tt_matrix(pz), tt_eye(pa.DOF_N, pa.occ))
    tt_tildeenergy = tt_kron(tt_tildeenergy, tt_tilenergy)
  
    #The total propogation Hamiltonian
    # Note that ficticious Harmonic oscillators carry negative sign
    H =  add_tensor(tt_He, tt_systemnumoc, small=eps) 
    H =  add_tensor(H, tt_tildenumoc, coeff = -1.0, small=eps) 
    H =  add_tensor(H, tt_systemenergy, coeff = 1.0, small=eps) 
    H =  add_tensor(H, tt_tildeenergy, coeff = 1.0, small=eps) 
    
    # Construct propagation operator, d/dt psi(t0)=-1j H psi(t0)
    H.nodes[0] *= -1j
    
    # convert to MPO
    A = MPS2MPO(H).truncation(small=eps)
    return A


def tt_eye(length,mode_dim):
    """
    generate Identity MPO (in the form of MPS) with rank-1.

    length: an INT number, specify the number of MPS nodes. 
    mode_dim: an INT number, specify the size of physical index.
    
    Returns: MPS
    -------
    """
    
    nb_arr = np.ones(length,dtype=int)*mode_dim
    y0 = MPS(length,nb_arr**2)
    
    for i in range(length):
        identy = np.eye(nb_arr[i],dtype=np.complex128)
        y0.nodes.append(identy.reshape((1,nb_arr[i]**2,1),order='F'))
    
    return y0


def tt_matrix(array):
    """
    generate a tensor train (in the form of MPS) from the input array.
    Now the array should be square matrix
    
    Returns: MPS
    -------
    """
    M, N = array.shape
    if(M!=N): 
        print("array shape not matched")
        sys.exit()
        
    y0 = MPS(1,np.array([M**2],dtype=int))
    y0.nodes.append(array.reshape((1,M**2,1),order='F'))
    
    return y0

def tt_kron(mps1,mps2):
    """
    kron product of two MPS

    Parameters
    ----------
    mps1 : TYPE MPS
    mps2 : TYPE MPS

    Returns
    -------
    New MPS which is the kron product of two MPS.

    """
    new_length = mps1.length + mps2.length
    
    new_nb = np.empty(new_length,dtype=int)
    
    new_mps = MPS(new_length,new_nb)
    
    for i in range(mps1.length):
        new_mps.nb[i] = mps1.nb[i]
        new_mps.nodes.append(mps1.nodes[i])
    for i in range(mps2.length):
        new_mps.nb[i+mps1.length] = mps2.nb[i]
        new_mps.nodes.append(mps2.nodes[i])
        
    return new_mps

def tt_ones(length,mode_dim):
    """
    generate MPS with rank-1, and all the elements are one.

    length: an INT number specify the number of MPS nodes. 
    mode_dim: an INT number specify the size of the physical index.
    
    Returns: MPS
    -------
    """
    
    nb_arr = np.ones(length,dtype=int)*mode_dim
    y0 = MPS(length,nb=nb_arr)  

    for i in range(length):
        tmp = np.ones((1,nb_arr[i],1),dtype=np.complex128)
        y0.nodes.append(tmp)

    return y0

def cal_property(mps0):
    """
    calculate the population and cohrerence for the spin-boson model
    input: mps0, the wavefunction MPS
    """
    sigma_arr = np.zeros(pa.DOF_E_SQ,dtype=np.complex128)
    
    su = pa.spin_up
    sd = pa.spin_down
    
    mps_up = mps0.copy()
    mps_down = mps0.copy()
    
    #ul = np.array([[1,0],[0,0]])
    ur = np.array([[0,1],[0,0]])
    mps_ur = mps0.copy()
    
    mps_up.nodes[0] = np.multiply(mps_up.nodes[0], su.reshape((1,-1,1)))
    mps_down.nodes[0] = np.multiply(mps_down.nodes[0], sd.reshape(1,-1,1))
    mps_ur.nodes[0] = (np.tensordot(ur, mps_ur.nodes[0],axes=((1),(1)))).transpose(1,0,2)
    sigma_arr[0] = calc_overlap(mps_up, mps_up)
    sigma_arr[3] = calc_overlap(mps_down, mps_down)
    sigma_arr[2] = calc_overlap(mps_up,mps_ur)
    sigma_arr[1] = calc_overlap(mps_ur,mps_up)
    return sigma_arr
    

def multiply_mps(mps1,mps2):
    """
    performing element wise multiplication of two MPS
    """
    
    nlen1 = len(mps1.nodes)
    nlen2 = len(mps2.nodes)
    if(nlen1!=nlen2):
        sys.exit('MPS length not matched! for multiply MPS')
    nb1 = mps1.nb
    nb2 = mps2.nb
    
    new_mps = MPS(nlen1,nb1)
    
    for i in range(nlen1):
        if(nb1[i]!=nb2[i]):
            sys.exit('MPS physical index not matched! for multiply MPS')
        ra1,ra2,ra3 = mps1.nodes[i].shape
        rb1,rb2,rb3 = mps2.nodes[i].shape
        
        array_tmp = np.einsum('ijk,ljm->iljkm',mps1.nodes[i],mps2.nodes[i]).reshape((ra1*rb1,nb1[i],ra3*rb3))
        
        new_mps.nodes.append(array_tmp)
    return new_mps.truncation(small=pa.eps)



def tt_tfd(initial_state, update_type='rk4', rk4slices = 1, mmax=4, RDO_arr_bench=None, show_steptime=False):

    """
    Performing the TT-TFD calculation
    
    update_type: In TDVP calculations, the method for updating each core of the MPS. 
                can be either 'rk4' or 'krylov'. can be 'rk4' or 'krylov'
    rk4slices: the timeslices for update_type='rk4'
    mmax: the size of Krylov subspace for update_type='krylov'
    RDO_arr_bench: An array. If given, print the error of each step's calculation compared to the benchmark RDO_arr_bench.
    show_steptime: if True, then show the time for each step. 
    
    """

    y0 = initial(initial_state)
    
    
    A = construct_Hamil(eps=pa.eps)

    RDO_arr = np.zeros((pa.TIME_STEPS, pa.DOF_E_SQ), dtype=np.complex128)
    t = np.arange(0, pa.TIME_STEPS * pa.DT, pa.DT)
 
    # Propagation loop
    START_TIME = time.time()
    print('Start doing propagation')
    
    for ii in range(pa.TIME_STEPS):
      
        print(ii,t[ii])
        STEP_TIME = time.time()
      
        #Doing TDVP
        y0 = tdvp1site(y0, A, pa.DT, update_type=update_type,mmax=mmax,rk4slices=rk4slices)
        #Calculating the Reduced density matrix
        RDO_arr[ii] = cal_property(y0)
      
        STEP_TIME2 = time.time()
        if(RDO_arr_bench is not None):
            compare_diff(RDO_arr[ii], RDO_arr_bench[ii]) 

        if(show_steptime):
            print('timefor tdvp:',STEP_TIME2-STEP_TIME)
      
    print("\tPropagation time:", time.time() - START_TIME)
  
    return t,RDO_arr 



def compare_diff(vec1,vec2):
    """
    testing function for comparing the two vectors
    """
    
    sr = 0.0; si = 0.0
    sr = np.max(np.abs(vec1.real-vec2.real))
    si = np.max(np.abs(vec1.imag-vec2.imag))
    print('error, real', sr, 'imag', si)
    #print('vec1',vec1)
    #print('vec2',vec2)

