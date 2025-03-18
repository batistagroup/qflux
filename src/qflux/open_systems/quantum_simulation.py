import numpy as np
import scipy.linalg as LA

#qiskit for quantum simulation
from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit.primitives import Estimator
from qiskit.quantum_info import SparsePauliOp

#the module in the current package
from . import trans_basis as tb
from . import dilation_circuit as dc
from .numerical_methods import Dynamics

#other utility functions for quantum simulation 
def expand(Gmat_org,Norg,Nexpand):
    """
    expand the propagator in the vectorized density matrix representation (e.g. match the 2^N dimension)
    """
    Gnew = np.zeros((Nexpand**2,Nexpand**2),dtype=np.complex128)
    for i in range(Norg):
        for j in range(Norg):
            for k in range(Norg):
                for l in range(Norg):
                    Gnew[i*Nexpand+j,k*Nexpand+l] = Gmat_org[i*Norg+j,k*Norg+l]
    return Gnew

def gen_Kraus_list(Gmat,N,tol=1E-5):
    """
      Generate the Kraus operators from the propagator with a given tolerance
      Input:
      - Gmat: matrix of the propagator (numpy array of shape (N^2, N^2)).
      - N: The system Hilbert space dimension
      - tol: tolerance for the Kraus operator representation.
      Returns:
      - Kraus: List of Kraus operators
    """
    #defining the Choi matrix from the matrix of the propagator
    C_mat = np.zeros(Gmat.shape,dtype=np.complex128)
    for i in range(N):
        for j in range(N):
            C_matij = np.zeros(Gmat.shape,dtype=np.complex128)
            for k in range(N):
                for l in range(N):
                    C_matij[i*N+k,l*N+j] = Gmat[j*N+k,l*N+i]
            C_mat += C_matij

    Kraus = []
    val,arr = LA.eigh(C_mat)
    for i in range(len(val)):
        if(val[i]>tol):
            Mi = np.sqrt(val[i])*arr[:,i].reshape(N,N)
            Kraus.append(Mi.conj().T)
    return Kraus

#class for quantum simulation
class DynamicsQ(Dynamics):
    

    def __init__(self, rep='Density', **kwargs):
        super().__init__(**kwargs)
        
        if(rep == 'Density'):
        #vectorized density matrix representation 
            self.rep    = 'Density'
            self.Nqb        = int(np.log2(self.Nsys**2))
        elif(rep == 'Kraus'):
        #Kraus operator representation 
            self.rep    = 'Kraus'
            self.Nqb        = int(np.log2(self.Nsys))
    
        #the counting qubits
        self.count_str = None
        #the observable 
        self.observable = None
        
        #the dilation methods for quantum simulation 
        self.dilation_method = 'Sz-Nagy'

    def set_dilation_method(self,method):
        """
        method: the dilation method, can be 'Sz-Nagy' or 'SVD' or 'SVD-Walsh'
        """
        self.dilation_method = method
        
    def set_count_str(self, count_str):
        self.count_str = count_str
            
    def set_observable(self, observable):
        self.observable = observable

    def init_statevec_vecdens(self):
        """
        get the dilated state vector from the initial density operator
        
        Returns: array of statevector (norm = 1), and initial norm
        -------
        """
        
        vec_rho0 = self.rho0.reshape(self.Nsys**2)
        norm0 = LA.norm(vec_rho0,2)
        
        #state vector in the dilated space
        statevec = vec_rho0/norm0
        
        return statevec, norm0

    def init_statevec_Kraus(self,tol=1E-6):
        """
        get the dilated state vector from the initial density operator 
            rho=\sum_i Pi |psi_i><psi_i|
        
        tol: ignore the state that have Pi<tol
        
        Returns: array of statevector |psi_i>, and probility Pi
        -------
        """
        
        #eigen decomposition of rho0 (density matrix is Hermitian)
        eigval,eigvec = LA.eigh(self.rho0)
        
        statevec = []
        prob = []
        for i in range(len(eigval)-1,-1,-1):
            if(abs(eigval[i])<tol): break
            prob.append(eigval[i])
            statevec.append(eigvec[:,i])
        
        return statevec, prob
    
    def _get_qiskit_observable(self, Isdilate = False, tol = 5E-3):
        
        if(self.observable is None): print('error, observable is None')
        
        if(Isdilate):
            num_qubits = self.Nqb+1 
            Obs_mat = np.zeros((2*self.Nsys,2*self.Nsys),dtype=np.complex128)
            Obs_mat[:self.Nsys,:self.Nsys] = self.observable[:self.Nsys,:self.Nsys]
        else:
            num_qubits = self.Nqb
            Obs_mat = self.observable

        Obs_paulis_dic = tb.ham_to_pauli(Obs_mat, num_qubits, tol=tol)
        
        #Prepare the qiskit observable from the pauli strings of observable matrix
        data = []
        coef = []
        for key in Obs_paulis_dic:
          data.append(key)
          coef.append(Obs_paulis_dic[key])
        obs_q = SparsePauliOp(data,coef)

        return obs_q 
    
    def qc_simulation_kraus(self, time_arr, shots=1024, Kraus = None, Gprop=None, tolk = 1E-5, tolo = 1E-5, **kwargs):
        """

        Parameters
        ----------
        time_arr : time array for simulation 
        shots : number of shots. The default is 1024.
        backend : Quantum Simulation backend. The default is AerSimulator().
        Kraus : Kraus operator dictionary (key: time_arr index, data: list of ndarray)
        Gprop : The propagator (may specified)
        tolk   : the tolerance of Kraus operator
        tolo   : the tolerance of decompose the observable to pauli matrix
        **kwargs: keyword arguments associate with the propagator calculation 

        Returns : array contain the quantum simulation result
        -------

        """
        nsteps = len(time_arr)
        
        #getting the Kraus operators 
        if(Kraus is None):
            Kraus = {}
            if(Gprop is None):
                print('calculating the propagator')
                Gprop = self.Gt_matrix_expo(time_arr,**kwargs)
            print('getting the Kraus operator')
            for i in range(nsteps):
                print('at i',i,'in total steps',nsteps)
                Kraus[i] = gen_Kraus_list(Gprop[i],self.Nsys,tol=tolk)
        print('calculating Kraus operator done')
        
        #performing the qiskit calculation 
        #Aer implementation of an Estimator
        estimator = Estimator()
        
        statevec, prob = self.init_statevec_Kraus()
        n_inistate = len(statevec)
        print('number of initial states in the density matrix', n_inistate)
        print('The probability',prob)
        
        #observable 
        obs_q = self._get_qiskit_observable(Isdilate = True, tol = tolo)
        
        print('starting the quantum simulation')
        result_simulation = np.zeros((nsteps),dtype=np.float64)
        
        for i in range(nsteps):
            print('istep',i,'in total', nsteps, 'steps')

            current_kraus_list = Kraus[i]
            print('number of Kraus opeartors',len(current_kraus_list))

            for ikraus in range(len(current_kraus_list)):
                
                for istate in range(n_inistate):
                    qc = self._create_circuit(current_kraus_list[ikraus], statevec[istate],Isscale=False)
                
                    result = estimator.run(qc, obs_q, shots = shots).result()
                    
                    result_simulation[i] += result.values[0]*prob[istate]
                    #print(result_simulation[i])
                    
        return result_simulation
        
        
    def qc_simulation_vecdens(self,time_arr,shots=1024, backend=AerSimulator(), Gprop = None, **kwargs):
        """

        Parameters
        ----------
        time_arr : time array for simulation 
        shots : number of shots. The default is 1024.
        backend : Quantum Simulation backend. The default is AerSimulator().
        Gprop : The propagator (may specified, or calculated when running)
        **kwargs: keyword arguments associate with the propagator calculation 

        Returns : array contain the quantum simulation result
        -------

        """
        
        if(Gprop is None):
            Gprop = self.Gt_matrix_expo(time_arr,**kwargs)
        
        nsteps = len(time_arr)
        
        if(self.count_str is None):  print("error: self.count_str unassigned")
        
        n_bitstr = len(self.count_str)
        
        statevec,norm0 = self.init_statevec_vecdens()
        
        result = np.zeros((nsteps,n_bitstr),dtype=np.float64)
        
        for i in range(nsteps):
            
            if(i%100==0):print('quantum simulation istep',i)
            
            Gt = Gprop[i]
            
            circuit,norm = self._create_circuit(Gt, statevec,Isscale=True)
            circuit.measure(range(self.Nqb+1),range(self.Nqb+1))
            if(self.dilation_method=='SVD-Walsh'):                
                circuit = transpile(circuit,backend)
           
            
            counts = backend.run(circuit,shots=shots).result().get_counts()
            
            for j in range(n_bitstr):
                bitstr = self.count_str[j]
                if(bitstr in counts):
                    result[i,j] = np.sqrt(counts[bitstr]/shots)*norm*norm0
                else:
                    print('at time', i,'shots=',shots,"no counts for", bitstr)
                
        return result
            
            
    def _create_circuit(self,array,statevec,Isscale=True):
        
        return dc.construct_circuit(self.Nqb, array, statevec, method = self.dilation_method, Isscale = Isscale)


