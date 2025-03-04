# Class for classical propagation methods
import numpy as np
import scipy.linalg as LA
import scipy.fft as sfft
from qutip import mesolve, Qobj

from . import params as pa
from . import trans_basis as tb


#Dynamics_grid
class Dynamics:
    '''
        Class for open-system dynamics (Lindblad equation). 

        Parameters
        ----------
        Nsys : int
            System Hilbert Space Dimension
        Hsys: 
            Hamiltonian of the system (numpy array of shape (N, N)).
        rho0: 
            Initial density matrix (numpy array of shape (N, N)).
        c_ops: 
            List of Collapse operators (numpy array of shape (N, N)).


    '''

    def __init__(self, Nsys, Hsys, rho0, c_ops = []):
        #--------- Required Attributes Populated During Execution ----------#
        self.Nsys           =    Nsys
        self.Hsys           =    Hsys
        self.rho0           =    rho0
        self.c_ops          =    c_ops


    def Gt_matrix_expo(self, time_arr, Is_show_step=False):
        """
            Getting the propagator of the Lindblad equation by matrix exponential
            Parameters:
            - H: Hamiltonian of the system (numpy array of shape (N, N)).
            - time_arr: Time array for dynamic simulation (array).
            - L: List of Collapse operators, with each operator is a numpy array of shape (N, N).
            Returns:
            - G_prop: List of propagators.
            """
        ident_h = np.eye(self.Nsys, dtype=np.complex128)
    
        # Amatrix for time-derivation of the vectorized density matrix
        Amat = -1j * (np.kron(self.Hsys, ident_h) - np.kron(ident_h, self.Hsys.T))
        for i in range(len(self.c_ops)):
          Amat += 0.5 * (2.0 * (np.kron(self.c_ops[i], self.c_ops[i].conj()))
                               - np.kron(ident_h, self.c_ops[i].T @ self.c_ops[i].conj())
                               - np.kron(self.c_ops[i].T.conj() @ self.c_ops[i], ident_h))
        
        G_prop = []
        for i in range(len(time_arr)):
            if(Is_show_step): print('step',i,'time',time_arr[i])
            Gt = LA.expm(Amat * time_arr[i])
            G_prop.append(Gt)
        return G_prop
    

    def propagate_matrix_exp(self, time_arr, observable, Is_store_state = False, Is_show_step = False, Is_Gt = False):
        """
            Solving the Lindblad equation by matrix exponential
            Parameters:
            - H: Hamiltonian of the system (numpy array of shape (N, N)).
            - rho0: Initial density matrix (numpy array of shape (N, N)).
            - time_arr: Time array for dynamic simulation (array).
            - L: List of Collapse operators, with each operator is a numpy array of shape (N, N).
            - observable: Observable for which the expectation value is computed (numpy array of shape (N, N)).
            - Is_store_state: Boolean variable that determines whether to output the density matrix list
            - show_step: Boolean variable that determines whether to print the current step during simulation
            Returns:
            - result: A class containing all the results
              result.expect: List of expectation values of the observable over time.
              result.G_prop: List of propagators.
              result.density_matrix: List of density matrices.
            """
            
        class Result:
          def __init__(self):
            self.expect = []
            if(Is_store_state):
              self.density_matrix = []
            if(Is_Gt):  self.Gprop = None
        result = Result()
    
        # Getting the propagator of the Lindblad equation 
        G_prop = self.Gt_matrix_expo(time_arr,Is_show_step)
        if(Is_Gt):  result.Gprop = G_prop
     
        # initialized vectorized density matrix
        vec_rho0 = self.rho0.reshape(self.Nsys**2)
            
        for i in range(len(time_arr)):
            
            vec_rhot = G_prop[i] @ vec_rho0
            
            # get the density matrix by reshaping
            rhot = vec_rhot.reshape(self.Nsys, self.Nsys)
    
            if(Is_store_state):  result.density_matrix.append(rhot)

            result.expect.append(np.trace(rhot @ observable).real)
    
        return result

            
    def propagate_qt(self, time_arr, observable, **kwargs):
        """
        Function used to propagate with qutip.mesolve function.
        
        - H: Hamiltonian of the system (Qobj).
        - rho0: Initial density matrix (Qobj).
        - time_arr: Time array for dynamic simulation (array).
        - c_ops: List of collapse operators (list of Qobj), can be empty for Liouville equation.
        - observable: Operators for which the expectation value is to be calculated .
        Returns:
        - expec_vals: List of expectation values of the observable over time.
        """
        c_ops = []
        for i in range(len(self.c_ops)):
            c_ops.append(Qobj(self.c_ops[i]))
        
        obs = []
        if(type(observable)==list):    
            for i in range(len(observable)):
                obs.append(Qobj(observable[i]))
        else:
            obs = Qobj(observable)
        
        result = mesolve(Qobj(self.Hsys), Qobj(self.rho0), time_arr, c_ops, obs, **kwargs)
        return result.expect


class DVR_grid:
    '''
        Class for Discrete Variable Representation (The system with the potential expressed in Grid point).  

        Parameters
        ----------
        Ngrid : int
            The number of grid points
        Hsys: 
            Hamiltonian of the system (numpy array of shape (N, N)).
        rho0: 
            Initial density matrix (numpy array of shape (N, N)).
        c_ops: 
            List of Collapse operators (numpy array of shape (N, N)).


    '''
   
    def __init__(self, xmin, xmax, Ngrid, mass):
        #--------- Required Attributes Populated During Execution ----------#
        self.Ngrid           =    Ngrid
        self.xmin          =    xmin
        self.xmax          =    xmax
        self.mass           =    mass 
        
        #associate with xgrid
        self.xgrid          =    None
        self._set_xgrid()
        self.dx             =    self.xgrid[1]-self.xgrid[0]
        
        #associate with kgrid
        self.dk             =    2.0*np.pi/((self.Ngrid)*self.dx)
        self.kgrid          =    None
        self.ak2            =    None #kinet energy array
        self.hamk           =    None #kinet hamiltonian
        self._set_kinet_ham()
               
        #potential
        self.potential      =    None 
        
        
        
        # self.rho0           =    rho0
        # self.c_ops          =    c_ops 
        
    def _set_xgrid(self):
        #set up the grid point
        self.xgrid = np.linspace(self.xmin,self.xmax,self.Ngrid)
        
    def set_potential(self,func_pot):
        """

        Setup the potential array
        ----------
        func_pot : a function f(x) that return the potential value at the point x

        """
        self.potential = np.zeros_like(self.xgrid)
        for i in range(self.Ngrid):
            self.potential[i] = func_pot(self.xgrid[i])
    
    def _set_kinet_ham(self):
        """
        setup kinetic hamiltonian Matrix in position x grid space
        """
        
        self.kgrid = np.zeros(self.Ngrid,dtype=np.float64)
        
        #ak2: kinetic energy array in k-space
        self.ak2   = np.zeros(self.Ngrid,dtype=np.float64)
        

        coef_k = pa.hbar**2/(2.0*self.mass)
        
        for i in range(self.Ngrid):
          if(i<self.Ngrid//2):
            self.kgrid[i] = i*self.dk
          else:
            self.kgrid[i] = -(self.Ngrid-i) * self.dk
        
          self.ak2[i] = coef_k*self.kgrid[i]**2

        akx0 = sfft.ifft(self.ak2)
        #hamk: kinetic hamiltonian Matrix in position x grid space
        self.hamk = np.zeros((self.Ngrid,self.Ngrid),dtype=np.complex128)
        
        for i in range(self.Ngrid):
          for j in range(self.Ngrid):
            if(i<j):
              self.hamk[i,j] = akx0[i-j].conj()
            else:
              self.hamk[i,j] = akx0[i-j]

        
    def get_eig_state(self,Nstate):
        """
        get the eigen state for potential in x-space
        Nstate: number of output eigen-state
        """
        Mata = self.hamk.copy()
        for i in range(self.Ngrid):
          Mata[i,i]+=self.potential[i]
      
        val,arr = LA.eigh(Mata)

        return val[:Nstate],arr[:,:Nstate]/self.dx**0.5
    
    def x2k_wave(self, psi):
        #transform the wavefunction from x space to k space
        psik = tb.x2k_wave(self.dx,psi)
        return psik
    
    def k2x_wave(self, psik):
        #transform the wavefunction from k space to x space
        psix  = tb.k2x_wave(self.dx,psik)
        return psix

        