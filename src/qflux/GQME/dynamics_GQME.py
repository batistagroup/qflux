"""
Class for GQME calculations 
"""

import numpy as np
import scipy.linalg as LA
from . import params as pa
from . import tt_tfd as tfd

import time
import sys


class DynamicsGQME:
    """
    Class for GQME calculation. 

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

    """
    def __init__(self, Nsys, Hsys, rho0):    
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

        
    
    def get_Liouvillian(self):
        """Get The projected Liouvillian from Hsys"""
        Isys = np.eye(self.Nsys)
        self.Liouv = np.kron(self.Hsys,Isys) -  np.kron(Isys,self.Hsys.T)
        
    def setup_propagator(self,Gt):
        """Setup the propagator from the input"""
        self.Gt = Gt        
    
    def setup_timestep(self,DT,TIME_STEPS):
        self.DT = DT
        self.TIME_STEPS = TIME_STEPS
        self.time_array = np.linspace(0,(TIME_STEPS-1)*DT,TIME_STEPS)

    def tt_tfd(self, initial_state=0, update_type='rk4', rk4slices = 1, mmax=4, RDO_arr_bench=None, show_steptime=False):
        """
        Performing the TT-TFD calculation
        
        update_type: In TDVP calculations, the method for updating each core of the MPS. 
                    can be either 'rk4' or 'krylov'. can be 'rk4' or 'krylov'
        rk4slices: the timeslices for update_type='rk4'
        mmax: the size of Krylov subspace for update_type='krylov'
        RDO_arr_bench: An array. If given, print the error of each step's calculation compared to the benchmark RDO_arr_bench.
        show_steptime: if True, then show the time for each step. 
        
        """
        time_arr, RDO_arr = tfd.tt_tfd(initial_state, update_type, rk4slices = rk4slices, mmax=mmax, RDO_arr_bench=RDO_arr_bench, show_steptime=show_steptime)
        return time_arr, RDO_arr
    
    def cal_propagator_tttfd(self):
        """
        Calculate the numerical exact propagator through TT-TFD method
        """
        
        U = np.zeros((pa.TIME_STEPS, pa.DOF_E_SQ, pa.DOF_E_SQ), dtype=np.complex128)
    
        # tt-tfd with initial state 0,1,2,3
        # initial state |0> means donor state |D>, |3> means acceptor state |A>
        # |1> is (|D> + |A>)/sqrt(2), |2> is (|D> + i|A>)/sqrt(2)
        print('========calculate the propagator, starting from 0 state========')
        t,U[:,:,0] = tfd.tt_tfd(0)
        print('========calculate the propagator, starting from 1 state========')
        t,U[:,:,1] = tfd.tt_tfd(1)
        print('========calculate the propagator, starting from 2 state========')
        t,U[:,:,2] = tfd.tt_tfd(2)
        print('========calculate the propagator, starting from 3 state========')
        t,U[:,:,3] = tfd.tt_tfd(3)
        print('========calculate the propagator done========')
    
        U_final = U.copy()
    
        # the coherence elements that start at initial state |D><A| and |A><D|
        # is the linear combination of above U results
        # |D><A| = |1><1| + i * |2><2| - 1/2 * (1 + i) * (|0><0| + |3><3|)
        U_final[:,:,1] = U[:,:,1] + 1.j * U[:,:,2] - 0.5 * (1. + 1.j) * (U[:,:,0] + U[:,:,3])
    
        # |A><D| = |1><1| - i * |2><2| - 1/2 * (1 - i) * (|0><0| + |3><3|)
        U_final[:,:,2] = U[:,:,1] - 1.j * U[:,:,2] - 0.5 * (1. - 1.j) * (U[:,:,0] + U[:,:,3])

        self.setup_propagator(U_final)
        
        return t,U_final
    
    def prop_puresystem(self):
        """Propogate the pure system"""
        Nstep = len(self.time_array)
        vec_rho = np.zeros((Nstep, self.Nsys**2), dtype=np.complex128)
        vec_rho[0] = self.vec_rho0.copy()
        for i in range(1,Nstep):
            vec_rho[i] = LA.expm(-1j*self.Liouv*(self.time_array[i]-self.time_array[i-1]))@vec_rho[i-1]
        return vec_rho
        
    def cal_F(self):
        """
        get the time-derivative of the propagator
        
        Returns
        -------
        F, Fdot: the 1st and 2nd order time-derivative of the propagator
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
    
        
    def _CalculateIntegral(self, F, linearTerm, prevKernel, kernel):
        """
        Function to Calculate Volterra Integral via Trapezoidal Rule
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
    
    def get_memory_kernel(self):
        """
        calculating the Memory kernel using volterra scheme
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


    def PropagateRK4(self,currentTime, memTime, kernel,
                 sigma_hold, sigma):

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
    
    def _Calculatef(self,currentTime, memTime, kernel, sigma_array, kVec):
    
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


    def solve_gqme(self, kernel, MEM_TIME, dtype = "Density"):
        """
        solve the GQME through RK4
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
