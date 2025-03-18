import numpy as np
from itertools import combinations, product
from scipy.linalg import expm, kron
from numpy.random import rand
random_numbers = np.random.uniform(low=10000, high=99999)
np.random.seed(int(random_numbers))
from .ansatz import *
from .ansatzVect import *

# Helper Functions for evolution
# def tag(A):
#     return A.tag if A else None #!!!!!!!a.any or a.all??
def tag(A):
    return getattr(A, 'tag', None) if A is not None else None


def aexp(A, theta):
    return expm(-1j * theta * A.mat / 2)

def add_A(A, P):
    A.A.append(P)
    A.theta = np.append(A.theta, 0)

def lmul(A, ψ):
    # print("ehrerw:",A)
    # print("amat:",np.shape(A.mat))
    # print("her in lmul",A.mat @ ψ)
    return A.mat @ ψ

def update_theta(ansatz, dtheta, dt):
    ansatz.theta = ansatz.theta + dtheta * dt

def partial_theta(A):
    len_A = len(A.A)
    ψ = A.ref.copy()
    res = []
    for i in range(len_A):
        ψ = aexp(A.A[i], A.theta[i]) @ ψ
        ψt = -0.5j * lmul(A.A[i], ψ)
        for j in range(i + 1, len_A):
            ψt = aexp(A.A[j], A.theta[j]) @ ψt
        res.append(ψt)
    return res

def build_m(ψ, dψ):
    l = len(dψ)
    res = np.zeros((l, l))
    ψ = ψ[:, np.newaxis] if len(ψ.shape) == 1 else ψ

    for μ in range(l):
        for ν in range(l):
            res[μ, ν] = 2 * np.real(dψ[μ].T.conj() @ dψ[ν] + ψ.T.conj() @ dψ[μ] @ ψ.T.conj() @ dψ[ν])
    return res

def update_m(m, dψa, ψ, dψ):
    l = m.shape[0]
    mp = np.zeros((l + 1, l + 1))
    dψa = dψa[:, np.newaxis]
    ψ = ψ[:, np.newaxis]
    mp[l, l] = 2 * np.real((dψa.T.conj() @ dψa + ψ.T.conj() @ dψa @ ψ.T.conj() @ dψa).item())
    mp[:l, :l] = m
    for μ in range(l):
        temp = 2 * np.real(dψ[μ].T.conj() @ dψa + ψ.T.conj() @ dψ[μ] @ ψ.T.conj() @ dψa)
        mp[μ, l] = temp.item()
        mp[l, μ] = temp.item()
    return mp

def build_v(ψ, dψ, He, Ha):
    L = -1j * (He - 1j * Ha)
    l = len(dψ)
    ψ = ψ[:, np.newaxis] if len(ψ.shape) == 1 else ψ
    res = np.zeros(l)
    for μ in range(l):
        res[μ] = 2 * np.real((dψ[μ].T.conj() @ L @ ψ + dψ[μ].T.conj() @ ψ @ ψ.T.conj() @ L.T.conj() @ ψ).item())
    return res

def update_v(v, dψₐ, ψ, He, Ha):
    L = -1j * (He - 1j * Ha)
    l = len(v)
    dψₐ = dψₐ[:, np.newaxis]
    ψ = ψ[:, np.newaxis]
    vp = np.zeros(l + 1)
    vp[:l] = v
    vp[l] = 2 * np.real((dψₐ.T.conj() @ L @ ψ + dψₐ.T.conj() @ ψ @ ψ.T.conj() @ L.T.conj() @ ψ).item())
    return vp

def lin_solve(m, v, λ=1e-4):
    A = m.T @ m
    try:
        xλ = np.linalg.solve(A + λ**2 * np.eye(A.shape[0]), m.T @ v)
    except np.linalg.LinAlgError:
        xλ = np.linalg.lstsq(A + λ**2 * np.eye(A.shape[0]), m.T @ v, rcond=None)[0]
    return xλ, 2 * v.T @ xλ - xλ.T @ m @ xλ

def get_newest_A(A):
    return A.A[-1] if A.A else None

def update_state(A):
    if A.A:
        state = A.ref.copy()
        for i in range(len(A.A)):
            state = aexp(A.A[i], A.theta[i]) @ state
        A.state = state
    else:
        A.state = A.ref

def set_ref(A, ref):
    A.ref = ref

def reset(A):
    A.theta = np.array([])
    A.A = []
    update_state(A)

class Result:
    def __init__(self, t, psi, energy, pop, theta, A, jump_t, jump_L, status):
        self.t = t
        self.psi = psi
        self.energy = energy  # Add energy attribute
        self.pop = pop
        self.theta = theta
        self.A = A
        self.jump_t = jump_t
        self.jump_L = jump_L
        self.status = status

    def __repr__(self):
        return f"Result(t={self.t}, psi={self.psi}, energy={self.energy}, pop={self.pop}, theta={self.theta}, A={self.A}, jump_t={self.jump_t}, jump_L={self.jump_L}, status={self.status})"

class AVQDSol:
    def __init__(self, t, u, θ, A, jump_L, jump_t, status=0, norm=[]):
        self.t = t
        self.u = u
        self.θ = θ
        self.A = A
        self.jump_L = jump_L
        self.jump_t = jump_t
        self.status = status
        self.norm = norm

def one_step(A, He, Ha, dt):

    relrcut = A.relrcut
    psi_ = A.state
    dpsi_ = partial_theta(A)
    M = build_m(psi_, dpsi_)
    V = build_v(psi_, dpsi_, He, Ha)
    dtheta, vmv = lin_solve(M, V)
    vmvMax = vmv
    Mtmp = M
    Vtmp = V
    dthetatmp = dtheta
    opTmp = None
    dpsi_ₐTmp = None
    add_flag = True

    while add_flag:
        # print("here in onestep While")
        tagTmp = None
        # print("Ansatz pool:",A.pool)
        for op in A.pool:
            # print("opertaro from pool:",(op.tag))
            if tag(get_newest_A(A)) == tag(op):
                # print("here in get newest_A")
                continue
            # print("here in op for loop")
            dpsi_ₐ = -0.5j * lmul(op, psi_)

            Mop = update_m(M, dpsi_ₐ, psi_, dpsi_)
            # print("Mop:",Mop)
            Vop = update_v(V, dpsi_ₐ, psi_, He, Ha)
            dthetaop, vmvOp = lin_solve(Mop, Vop)
            # print("berfore the if:","vmvOp:",vmvOp, "vmvMax:",vmvMax)

            if vmvOp > vmvMax:
                # print("vmvOp:",vmvOp, "vmvMax:",vmvMax)

                Mtmp = Mop
                Vtmp = Vop
                vmvMax = vmvOp
                dthetatmp = dthetaop
                opTmp = op
                tagTmp = tag(op)
                # print("tagTmp:",tagTmp)
                dpsi_ₐTmp = dpsi_ₐ
        # print("vmvMax - vmv:",vmvMax-vmv)
        add_flag = vmvMax - vmv >= relrcut

        

        if tagTmp is not None and add_flag:

            add_A(A, opTmp)
            vmv = vmvMax
            M = Mtmp
            V = Vtmp
            dtheta = dthetatmp
            dpsi_.append(dpsi_ₐTmp)

    update_theta(A, dtheta, dt)
    # print("Ansatz:",A.theta)
    update_state(A)


def solve_avq_traj(H, A, tspan, dt, save_state=True, save_everystep=True):
    # Store initial reference state for reinitialization
    ref_init = A.ref.copy()
    update_state(A)
    
    # Initialize time and recorder lists
    t = tspan[0]
    t_list = []
    u_list = []
    theta_list = []
    A_list = []
    jump_L_list = []
    jump_t_list = []

    # Break flag for error handling
    break_flag = False
    Gamma = 0  # Jump recorder
    q = rand()  # Random number for quantum jump
    # print("expo:",np.exp(-Gamma))
    # print("random q:", q)

    if save_everystep:
        t_list.append(t)
        theta_list.append(A.theta.copy())

        A_list.append([tag(a) for a in A.A])
        if save_state:
            u_list.append(A.state.copy())

    # Main evolution loop
    # This loop performs the simulation of quantum system evolution over a set time span.
    # It integrates the system state forward in time while handling both deterministic 
    # evolution and stochastic quantum jumps.

    while t + dt <= tspan[1]:  # Loop over the total time span with step size dt
        # print("timestep:", dt)  # Print current time step size for monitoring
        He = H.He  # Extract Hermitian Hamiltonian component for deterministic evolution
        Ha = H.Ha  # Extract Antihermitian Hamiltonian component for jump handling

        if np.exp(-Gamma) > q:  # Check if no quantum jump occurs, based on jump probability threshold
            try:
                one_step(A, He, Ha, dt)  # Perform deterministic evolution for one time step
            except Exception:  # In case of failure, exit the loop with a flag
                break_flag = True
                break
            psi_ = A.state  # Update the current state of the system
            Gamma += 2 * np.real(psi_.T.conj() @ Ha @ psi_) * dt  # Update accumulated Gamma for jump probability

            if save_everystep:  # Option to save state at every time step
                t_list.append(t + dt)  # Save current time
                theta_list.append(A.theta.copy())  # Save current ansatz parameters
                A_list.append([tag(a) for a in A.A])  # Save current state of ansatz components

        else:  # Quantum jump occurs
            psi_ = A.state  # Get current state before jump

            # Compute the probabilities for different jump operators
            weights = np.array([np.real(psi_.T.conj() @ LdL @ psi_) for LdL in H.LdL])

            # Randomly choose a jump operator based on computed probabilities
            L_index = np.random.choice([i for i in range(np.shape(H.Llist)[0])], p=weights/weights.sum())
            L = H.Llist[L_index]  # Select corresponding jump operator

            # Apply the jump operator to the state and normalize it
            psi_n = L @ psi_ / np.linalg.norm(L @ psi_)

            # Record the jump event details (time and operator used)
            jump_t_list.append(t)
            jump_L_list.append(tag(L))

            # Record the state and parameters post-jump
            t_list.append(t + dt)
            set_ref(A, psi_n)  # Update the state with post-jump state
            reset(A)  # Reset the ansatz parameters after jump

            # Save updated ansatz parameters post-jump
            theta_list.append(A.theta.copy())
            A_list.append([tag(a) for a in A.A])

            # Reset the jump probability accumulator
            Gamma = 0
            q = rand()  # Generate a new random number for the next jump check

        # Optionally save the state at each time step
        if save_state:
            u_list.append(A.state.copy())  # Store the current state

        # Increment the simulation time by dt
        t += dt


    # Final check for status if no errors occurred
    if not break_flag:
        if (not jump_t_list or not np.isclose(t, jump_t_list[-1])) and not save_everystep:
            theta_list.append(A.theta.copy())
            A_list.append([tag(a) for a in A.A])

    # Reset the ansatz to the initial condition
    # print("ansatz before reset:",A.A)
    set_ref(A, ref_init)
    reset(A)
    # print("tlist:",t_list)
    # Return solution object
    return AVQDSol(t_list, u_list, theta_list, A_list, jump_L_list, jump_t_list, status=0 if not break_flag else 1)

# Function to evolve a single trajectory and collect results
def solve_avq_trajectory(H, ansatz, tf, dt):
    # Solve for a single trajectory using solve_avq
    res = solve_avq_traj(H, ansatz, [0, tf], dt, save_state=True, save_everystep=True)
    
    # Post-process the results to extract energy and populations
    tlist = res.t
    psi_list = res.u
    energy = []
    pop = []
    
    for i in range(len(tlist)):
        psi = psi_list[i]
        energy.append(np.real(np.conj(psi.T) @ (H.He) @ psi))  # Energy
        pop.append([np.abs(psi[0])**2, np.abs(psi[1])**2])  # Population in ground/excited states for amplitude damping

    return Result(tlist, psi_list, energy, pop, res.θ, res.A, res.jump_t, res.jump_L, res.status)
    pass

def solve_avq_vect(H, A, tspan, dt):
    ref_init = A.ref.copy()
    nqbit = A.nqbit
    update_state(A)
    t = tspan[0]
    t_list = [t]
    u_list = [A.state.reshape(2**nqbit, 2**nqbit)]
    theta_list = [A.theta.copy()]
    A_list = [[tag(a) for a in A.A]]
    norm_list = [1.0]
    Gamma = 0
    while t + dt <= tspan[1]:
        He = H.He
        Ha = H.Ha
        one_step(A, He, Ha, dt)
        psi_ = A.state
        Gamma += 2 * np.real(psi_.T.conj() @ Ha @ psi_) * dt
        ρ = psi_.reshape(2**nqbit, 2**nqbit)
        ρ /= np.trace(ρ)
        t += dt
        t_list.append(t)
        u_list.append(ρ)
        theta_list.append(A.theta.copy())
        A_list.append([tag(a) for a in A.A])
        norm_list.append(np.exp(-Gamma))
    set_ref(A, ref_init)
    reset(A)
    return AVQDSol(t_list, u_list, theta_list, A_list, [], [], norm=norm_list)
