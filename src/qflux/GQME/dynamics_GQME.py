"""
Class for GQME calculations 
"""

import numpy as np
import scipy.linalg as LA
from . import params as pa
from . import tt_tfd as tfd
from typing import Tuple, Optional

import time
import sys


class DynamicsGQME:
    """
    Class for Generalized Quantum Master Equation (GQME) calculation.
    
    This class provides methods for calculating the memory kernel, 
    performing TT-TFD calculations, and solving the GQME for open quantum systems.
    
    Attributes:
        Nsys (int): System Hilbert space dimension.
        Hsys (np.ndarray): System Hamiltonian of shape (N, N).
        rho0 (np.ndarray): Initial density matrix of shape (N, N).
        vec_rho0 (np.ndarray): Vectorized initial density matrix of shape (N^2,).
        Liouv (Optional[np.ndarray]): Liouvillian superoperator.
        DT (Optional[float]): Time step for evolution.
        TIME_STEPS (Optional[int]): Number of time steps.
        time_array (Optional[np.ndarray]): Array of time points.
        Gt (Optional[np.ndarray]): Time-evolution propagator.
    """
    
    def __init__(self, Nsys: int, Hsys: np.ndarray, rho0: np.ndarray) -> None:
        """
        Initialize the DynamicsGQME instance.
        
        Args:
            Nsys (int): System Hilbert space dimension.
            Hsys (np.ndarray): System Hamiltonian of shape (N, N).
            rho0 (np.ndarray): Initial density matrix of shape (N, N).
        """
        self.Nsys           =    Nsys
        self.Hsys           =    Hsys
        self.rho0           =    rho0
        self.vec_rho0       =    rho0.reshape(Nsys**2)
        
        #Liouvillian
        self.Liouv          =   None
        self.get_Liouvillian()    
        
        #time 
        self.DT             =   None
        self.TIME_STEPS     =   None
        self.time_array     =   None
        
        #propagator
        self.Gt             =   None

    
    def get_Liouvillian(self) -> None:
        r"""
        Construct the Liouvillian superoperator using the system Hamiltonian.
        
        The Liouvillian is defined as:
            L = H \otimes I - I \otimes H^T
        """
        Isys = np.eye(self.Nsys)
        self.Liouv = np.kron(self.Hsys,Isys) -  np.kron(Isys,self.Hsys.T)
        
    def setup_propagator(self, Gt: np.ndarray) -> None:
        """
        Set the time-evolution propagator.

        Args:
            Gt (np.ndarray): Time-dependent propagator.
        """
        self.Gt = Gt        
    
    def setup_timestep(self, DT: float, TIME_STEPS: int) -> None:
        """
        Set up the time discretization parameters.

        Args:
            DT (float): Time step size.
            TIME_STEPS (int): Number of time steps.
        """
        self.DT = DT
        self.TIME_STEPS = TIME_STEPS
        self.time_array = np.linspace(0,(TIME_STEPS-1)*DT,TIME_STEPS)

    def tt_tfd(
        self,
        initial_state: int = 0,
        update_type: str = 'rk4',
        rk4slices: int = 1,
        mmax: int = 4,
        RDO_arr_bench: Optional[np.ndarray] = None,
        show_steptime: bool = False
        ) -> Tuple[np.ndarray, np.ndarray]:
        
        """
        Perform Tensor-Train Thermofield Dynamics (TT-TFD) calculation.
        
        Args:
            initial_state (int, optional): Index of the initial state. Defaults to 0.
            update_type (str, optional): Method for time evolution. Either 'rk4' or 'krylov'.
            rk4slices (int, optional): Number of RK4 substeps if update_type is 'rk4'. Defaults to 1.
            mmax (int, optional): Dimension of Krylov subspace if update_type is 'krylov'. Defaults to 4.
            RDO_arr_bench (Optional[np.ndarray], optional): Benchmark RDO array to compute error at each step. Defaults to None.
            show_steptime (bool, optional): Whether to print timing for each step. Defaults to False.
        
        Returns:
            Tuple[np.ndarray, np.ndarray]: Time array and RDO (reduced density operator) array.
        """
        time_arr, RDO_arr = tfd.tt_tfd(
            initial_state,
            update_type,
            rk4slices=rk4slices,
            mmax=mmax,
            RDO_arr_bench=RDO_arr_bench,
            show_steptime=show_steptime
        )

        return time_arr, RDO_arr
    
    def cal_propagator_tttfd(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate the numerical exact propagator for Spin-Boson model using the TT-TFD method.

        Returns:
            Tuple[np.ndarray, np.ndarray]: Time array and 3D propagator array.
        """
        
        U = np.zeros((pa.TIME_STEPS, pa.DOF_E_SQ, pa.DOF_E_SQ), dtype=np.complex128)
    
        # tt-tfd with initial state 0,1,2,3
        # initial state |0> means donor state |D>, |3> means acceptor state |A>
        # |1> is (|D> + |A>)/sqrt(2), |2> is (|D> + i|A>)/sqrt(2)
        for i in range(4):
            print(f"======== calculate the propagator, starting from state {i} ========")
            t, U[:, :, i] = tfd.tt_tfd(i)
   
        U_final = U.copy()
    
        # the coherence elements that start at initial state |D><A| and |A><D|
        # is the linear combination of above U results
        # |D><A| = |1><1| + i * |2><2| - 1/2 * (1 + i) * (|0><0| + |3><3|)
        U_final[:,:,1] = U[:,:,1] + 1.j * U[:,:,2] - 0.5 * (1. + 1.j) * (U[:,:,0] + U[:,:,3])
    
        # |A><D| = |1><1| - i * |2><2| - 1/2 * (1 - i) * (|0><0| + |3><3|)
        U_final[:,:,2] = U[:,:,1] - 1.j * U[:,:,2] - 0.5 * (1. - 1.j) * (U[:,:,0] + U[:,:,3])

        self.setup_propagator(U_final)
        print('========calculate the propagator done========')
        
        return t,U_final
    
    def prop_puresystem(self) -> np.ndarray:
        """
        Propagate the pure system under the unitary evolution.

        Returns:
            np.ndarray: Vectorized density matrix at each time step (shape (TIME_STEPS, N^2)).
        """
        Nstep = len(self.time_array)
        vec_rho = np.zeros((Nstep, self.Nsys**2), dtype=np.complex128)
        vec_rho[0] = self.vec_rho0.copy()
        for i in range(1,Nstep):
            vec_rho[i] = LA.expm(-1j*self.Liouv*(self.time_array[i]-self.time_array[i-1]))@vec_rho[i-1]
        return vec_rho
        
    def cal_F(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute the first and second order time derivatives of the propagator.

        Returns:
            Tuple[np.ndarray, np.ndarray]: F and Fdot matrices (each of shape (TIME_STEPS, N^2, N^2)).
        """

        F = np.zeros((self.TIME_STEPS, self.Nsys**2, self.Nsys**2), dtype=np.complex128)
        Fdot = np.zeros((self.TIME_STEPS, self.Nsys**2, self.Nsys**2), dtype=np.complex128)
        
        if(self.Gt is None):
            print('error in cal_F: please setup the pre-calculated propagator super-operator')
            sys.exit()
            
        for j in range(self.Nsys**2):
            for k in range(self.Nsys**2):
                # extracts real and imag parts of U element
                Ureal = self.Gt[:,j,k].copy().real
                Uimag = self.Gt[:,j,k].copy().imag
    
                # F = i * d/dt U so Re[F] = -1 * d/dt Im[U] and Im[F] = d/dt Re[U]
                Freal = -1. * np.gradient(Uimag.flatten(), self.DT, edge_order = 2)
                Fimag = np.gradient(Ureal.flatten(), self.DT, edge_order = 2)
    
                # Fdot = d/dt F so Re[Fdot] = d/dt Re[F] and Im[Fdot] = d/dt Im[F]
                Fdotreal = np.gradient(Freal, self.DT)
                Fdotimag = np.gradient(Fimag, self.DT)
    
                F[:,j,k] = Freal[:] + 1.j * Fimag[:]
                Fdot[:,j,k] = Fdotreal[:] + 1.j * Fdotimag[:]
               
        return F,Fdot
    
        
    def _CalculateIntegral(
        self,
        F: np.ndarray,
        linearTerm: np.ndarray,
        prevKernel: np.ndarray,
        kernel: np.ndarray
        ) -> np.ndarray:
        """
        Compute the Volterra integral using the trapezoidal rule.

        Args:
            F (np.ndarray): Derivative of the propagator.
            linearTerm (np.ndarray): Linear term from memory kernel equation.
            prevKernel (np.ndarray): Kernel from previous iteration.
            kernel (np.ndarray): Kernel to be updated.

        Returns:
            np.ndarray: Updated kernel.
        """
        
        # time step loop starts at 1 because K is equal to linear part at t = 0
        for n in range(1, self.TIME_STEPS):
            kernel[n,:,:] = 0.
    
            # f(a) and f(b) terms
            kernel[n,:,:] += 0.5 * self.DT * F[n,:,:] @ kernel[0,:,:]
            kernel[n,:,:] += 0.5 * self.DT * F[0,:,:] @ prevKernel[n,:,:]
    
            # sum of f(a + kh) term
            for c in range(1, n):
                # since a new (supposed-to-be-better) guess for the
                # kernel has been calculated for previous time steps,
                # can use it rather than prevKernel
                kernel[n,:,:] += self.DT * F[n - c,:,:] @ kernel[c,:,:]
    
            # multiplies by i and adds the linear part
            kernel[n,:,:] = 1.j * kernel[n,:,:] + linearTerm[n,:,:]
    
        return kernel
    
    def get_memory_kernel(self) -> np.ndarray:
        """
        Compute the memory kernel using the Volterra scheme.

        Returns:
            np.ndarray: Memory kernel array (shape (TIME_STEPS, N^2, N^2)).
        """
        F,Fdot = self.cal_F()
        
        linearTerm = 1.j * Fdot.copy() # first term of the linear part
        for l in range(self.TIME_STEPS):
            # subtracts second term of linear part
            linearTerm[l,:,:] -= 1./pa.HBAR * F[l,:,:] @ self.Liouv
            
        START_TIME = time.time() # starts timing
        # sets initial guess to the linear part
        prevKernel = linearTerm.copy()
        kernel = linearTerm.copy()

        # loop for iterations
        for numIter in range(1, pa.MAX_ITERS + 1):
        
            iterStartTime = time.time() # starts timing of iteration
            print("Iteration:", numIter)
        
            # calculates kernel using prevKernel and trapezoidal rule
            kernel = self._CalculateIntegral(F, linearTerm, prevKernel, kernel)
        
            numConv = 0 # parameter used to check convergence of entire kernel
            for i in range(self.Nsys**2):
                for j in range(self.Nsys**2):
                    for n in range(self.TIME_STEPS):
                        # if matrix element and time step of kernel is converged, adds 1
                        if abs(kernel[n][i][j] - prevKernel[n][i][j]) <= pa.CONVERGENCE_PARAM:
                            numConv += 1
        
                        # if at max iters, prints which elements and time steps did not
                        # converge and prevKernel and kernel values
                        elif numIter == pa.MAX_ITERS:
                            print("\tK time step and matrix element that didn't converge: %s, %s%s"%(n,i,j))
        
            print("\tIteration time:", time.time() - iterStartTime)
        
            # enters if all times steps and matrix elements of kernel converged
            if numConv == self.TIME_STEPS * self.Nsys**2 * self.Nsys**2:
                # prints number of iterations and time necessary for convergence
                print("Number of Iterations:", numIter, "\tVolterra time:", time.time() - START_TIME)
        
                break # exits the iteration loop
        
            # if not converged, stores kernel as prevKernel, zeros the kernel, and then
            # sets kernel at t = 0 to linear part
            prevKernel = kernel.copy()
            kernel = linearTerm.copy()
        
            # if max iters reached, prints lack of convergence
            if numIter == pa.MAX_ITERS:
                print("\tERROR: Did not converge for %s iterations"%pa.MAX_ITERS)
                print("\tVolterra time:", print(time.time() - START_TIME))
        
        return kernel


    def PropagateRK4(
        self,
        currentTime: float,
        memTime: float,
        kernel: np.ndarray,
        sigma_hold: np.ndarray,
        sigma: np.ndarray
        ) -> np.ndarray:
        """
        Perform one 4th-order Runge-Kutta (RK4) integration step for GQME.
    
        Args:
            currentTime (float): Current time.
            memTime (float): Memory time cutoff.
            kernel (np.ndarray): Memory kernel array (shape (TIME_STEPS, N^2, N^2)).
            sigma_hold (np.ndarray): Current vectorized reduced density matrix (N^2,).
            sigma (np.ndarray): Array of vectorized reduced density matrix up to current time (TIME_STEPS, N^2).
    
        Returns:
            np.ndarray: Updated vectorized reduced density matrix after one RK4 step (N^2,).
        """
        f_0 = self._Calculatef(currentTime, memTime,
                         kernel, sigma, sigma_hold)
    
        k_1 = sigma_hold + self.DT * f_0 / 2.
        f_1 = self._Calculatef(currentTime + self.DT / 2., memTime,
                         kernel, sigma, k_1)
    
        k_2 = sigma_hold + self.DT * f_1 /2.
        f_2 = self._Calculatef(currentTime + self.DT / 2., memTime,
                         kernel, sigma, k_2)
    
        k_3 = sigma_hold + self.DT * f_2
        f_3 = self._Calculatef(currentTime + self.DT, memTime,
                         kernel, sigma, k_3)
    
        sigma_hold += self.DT / 6. * (f_0 + 2. * f_1 + 2. * f_2 + f_3)
    
        return sigma_hold
    
    def _Calculatef(
                    self,
                    currentTime: float,
                    memTime: float,
                    kernel: np.ndarray,
                    sigma_array: np.ndarray,
                    kVec: np.ndarray
                ) -> np.ndarray:
        """
        Evaluate the time derivative in GQME.
        
        Args:
            currentTime (float): Current time.
            memTime (float): Memory time cutoff.
            kernel (np.ndarray): Memory kernel (shape (TIME_STEPS, N^2, N^2)).
            sigma_array (np.ndarray): History of vectorized system density matrix (TIME_STEPS, N^2).
            kVec (np.ndarray): Vector to be evolved (N^2,).
        
        Returns:
            np.ndarray: Time derivative vector f(t) (N^2,).
        """
        memTimeSteps = int(memTime / self.DT)
        currentTimeStep = int(currentTime / self.DT)
    
        f_t = np.zeros(kVec.shape, dtype=np.complex128)
    
        f_t -= 1.j / pa.HBAR * self.Liouv @ kVec
    
        limit = memTimeSteps
        if currentTimeStep < (memTimeSteps - 1):
            limit = currentTimeStep
        for l in range(limit):
            f_t -= self.DT * kernel[l,:,:] @ sigma_array[currentTimeStep - l]
    
        return f_t


    def solve_gqme(
                self,
                kernel: np.ndarray,
                MEM_TIME: float,
                dtype: str = "Density"
            ) -> np.ndarray:
        """
        Solve the GQME using RK4 integration.
    
        Args:
            kernel (np.ndarray): Memory kernel (shape (TIME_STEPS, N^2, N^2)).
            MEM_TIME (float): Memory cutoff time.
            dtype (str): Type of data to propagate: "Density" or "Propagator".
    
        Returns:
            np.ndarray: Propagated state (shape depends on dtype).
        """
        
        if(dtype=="Density"):
            # array for reduced density matrix elements
            sigma = np.zeros((self.TIME_STEPS, self.Nsys**2), dtype=np.complex128)
            # array to hold copy of sigma
            sigma_hold = np.zeros(self.Nsys**2, dtype = np.complex128)
            
            # sets the initial state
            sigma[0,:] = self.vec_rho0.copy()
            sigma_hold = self.vec_rho0.copy()
        elif(dtype == "Propagator"):
            # array for reduced density matrix elements
            sigma = np.zeros((self.TIME_STEPS, self.Nsys**2, self.Nsys**2), dtype=np.complex128)
            # array to hold copy of sigma
            sigma_hold = np.zeros(self.Nsys**2, dtype = np.complex128)
            
            #time 0 propagator: identity superoperator
            sigma[0] = np.eye(self.Nsys**2)
            #array to hold copy of G propagator
            sigma_hold = np.eye((self.Nsys**2), dtype=np.complex128)
        else:
            sys.exit('GQME input error, dtype should be "Density" or "Propagator"')

        # loop to propagate sigma
        print(">>> Starting GQME propagation, memory time =", MEM_TIME)
        for l in range(self.TIME_STEPS - 1): # it propagates to the final time step
            if l%100==0: print(l)
            currentTime = l * self.DT
        
            sigma_hold = self.PropagateRK4(currentTime, MEM_TIME, kernel, sigma_hold, sigma)
        
            sigma[l + 1] = sigma_hold.copy()
        
        return sigma
