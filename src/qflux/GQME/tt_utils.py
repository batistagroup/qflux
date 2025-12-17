from .params_class import Params
import numpy as np 
from mpsqd.utils import MPS, add_tensor, MPS2MPO, calc_overlap
from .tt_tfd import (initial, construct_Hamil, tt_eye, tt_matrix, 
                     tt_kron, tt_ones, multiply_mps, compare_diff,
                     cal_property, tt_tfd)
from typing import Any, Callable
                        
pp = Params()

def tt_initial_state(istate: int, pp=pp) -> MPS:
    """
    Initialize the state in tensor-train (MPS) format for a TT-TFD calculation.

    Parameters
    ----------
    istate : int
        Type of initial electronic state:
        0 : spin-up
        1 : (spin-up + spin-down) / sqrt(2)
        2 : (spin-up + i * spin-down) / sqrt(2)
        3 : spin-down

    Returns
    -------
    MPS
        Initialized MPS with the chosen electronic state at the first site
        and vacuum/ground states on the remaining sites.
        QFlux uses mpsqd https://github.com/qiangshi-group/MPSQD
    """
    # Sanity check on istate
    if istate not in (0, 1, 2, 3):
        raise ValueError(f"Invalid istate={istate}. Must be in {{0, 1, 2, 3}}.")

    # -------------------------------------------------------------------------
    # Define single-site electronic tensors
    # -------------------------------------------------------------------------
    su = np.zeros((1, pp.DOF_E, pp.MAX_TT_RANK), dtype=np.complex128)
    sd = np.zeros((1, pp.DOF_E, pp.MAX_TT_RANK), dtype=np.complex128)

    su[0, :, 0] = pp.spin_up
    sd[0, :, 0] = pp.spin_down

    # Superpositions
    inv_sqrt2 = 1.0 / np.sqrt(2.0)
    e1 = inv_sqrt2 * (su + sd)
    e2 = inv_sqrt2 * (su + 1j * sd)

    # Select the initial electronic core
    electronic_cores = {
        0: su,
        1: e1,
        2: e2,
        3: sd,
    }
    first_core = electronic_cores[istate]

    # -------------------------------------------------------------------------
    # Build MPS structure
    # -------------------------------------------------------------------------
    # nbarr: local dimensions for each site
    num_sites = 1 + 2 * pp.DOF_N
    nbarr = np.full(num_sites, pp.occ, dtype=int)
    nbarr[0] = pp.DOF_E  # first site is electronic

    y0 = MPS(num_sites, nb=nbarr)
    y0.nodes.append(first_core)

    # Middle sites: identity-like / vacuum cores
    middle_core = np.zeros(
        (pp.MAX_TT_RANK, pp.occ, pp.MAX_TT_RANK),
        dtype=np.complex128
    )
    middle_core[0, 0, 0] = 1.0

    # Append 2 * DOF_N - 1 middle cores
    for _ in range(2 * pp.DOF_N - 1):
        y0.nodes.append(middle_core)

    # Last site: right boundary core with rank-1 right bond
    last_core = np.zeros(
        (pp.MAX_TT_RANK, pp.occ, 1),
        dtype=np.complex128
    )
    last_core[0, 0, 0] = 1.0
    y0.nodes.append(last_core)

    return y0


def discretize_ohmic(freq_count: int, pp=pp):
    """
    Discretize an Ohmic spectral density into `freq_count` modes.

    Returns
    -------
    freq : (N,) array
    ck   : (N,) array
    gk   : (N,) array
    thetak, sinhthetak, coshthetak : (N,) arrays
    """
    N = freq_count

    om = pp.OMEGA_C / N * (1.0 - np.exp(-pp.OMEGA_MAX / pp.OMEGA_C))

    freq = np.zeros(N, dtype=float)
    ck = np.zeros(N, dtype=float)
    gk = np.zeros(N, dtype=float)
    thetak = np.zeros(N, dtype=float)
    sinhthetak = np.zeros(N, dtype=float)
    coshthetak = np.zeros(N, dtype=float)

    for i in range(N):
        freq[i] = -pp.OMEGA_C * np.log(
            1.0 - (i + 1) * om / pp.OMEGA_C
        )
        ck[i] = np.sqrt(pp.XI * om) * freq[i]
        gk[i] = -ck[i] / np.sqrt(2.0 * freq[i])

        th = np.arctanh(np.exp(-pp.BETA * freq[i] / 2.0))
        thetak[i] = th
        sinhthetak[i] = np.sinh(th)
        coshthetak[i] = np.cosh(th)

    return freq, ck, gk, thetak, sinhthetak, coshthetak


def build_electronic_hamiltonian(epsilon: float, gamma_da: float):
    """2x2 electronic Hamiltonian in matrix form."""
    px = np.array([[0.0, 1.0],
                   [1.0, 0.0]], dtype=np.complex128)
    pz = np.array([[1.0, 0.0],
                   [0.0, -1.0]], dtype=np.complex128)
    return epsilon * pz + gamma_da * px


def tt_embed_electronic(tt_He, total_boson_modes: int, occ: int):
    """
    Kronecker-extend a 2x2 TT-matrix to include 2*DOF_N bosonic modes.
    """
    return tt_kron(tt_He, tt_eye(2 * total_boson_modes, occ))


def build_number_operator_local(occ: int):
    """Local harmonic number operator in matrix form."""
    return np.diag(np.arange(occ, dtype=np.complex128))


def tt_zero_like_eye(num_sites: int, occ: int):
    """Create a TT with identity structure and then zero all cores."""
    tt_obj = tt_eye(num_sites, occ)
    for i in range(num_sites):
        tt_obj.nodes[i] *= 0.0
    return tt_obj


def tt_sum_local_operators(num_sites: int,
                           occ: int,
                           local_mats,
                           site_coeffs,
                           eps: float):
    """
    Build sum_k site_coeffs[k] * (I ... x local_mats[k] x ... I) in TT form.
    """
    tt_total = tt_zero_like_eye(num_sites, occ)

    for k, (Mloc, coeff) in enumerate(zip(local_mats, site_coeffs)):
        tmp0 = tt_matrix(Mloc)
        tmp0.nodes[0] *= coeff

        if k == 0:
            tmp = tt_kron(tmp0, tt_eye(num_sites - 1, occ))
        elif k < num_sites - 1:
            tmp = tt_kron(tt_eye(k - 1, occ), tmp0)
            tmp = tt_kron(tmp, tt_eye(num_sites - k, occ))
        else:  # last site
            tmp = tt_kron(tt_eye(k, occ), tmp0)

        tt_total = add_tensor(tt_total, tmp, small=eps)

    return tt_total


def build_displacement_local(occ: int):
    """
    Local displacement operator (x operator) in HO basis.
    """
    D = np.zeros((occ, occ), dtype=np.complex128)
    for i in range(occ - 1):
        s = np.sqrt(i + 1.0)
        D[i, i + 1] = s
        D[i + 1, i] = s
    return D


def tt_number_operator_physical(freq, eps: float, pp=pp):
    r"""
    TT representation of sum_k freq[k] * a_k^\dagger a_k on DOF_N sites.
    """
    N = pp.DOF_N
    numoc = build_number_operator_local(pp.occ)
    local_mats = [numoc] * N
    return tt_sum_local_operators(N, pp.occ, local_mats, freq, eps)


def tt_displacement_physical(gk, coshthetak, eps: float, pp=pp):
    r"""
    TT representation of sum_k gk[k] cosh(theta_k) (a_k + a_k^\dagger)
    """
    N = pp.DOF_N
    D = build_displacement_local(pp.occ)
    local_mats = [D] * N
    coeffs = gk * coshthetak
    return tt_sum_local_operators(N, pp.occ, local_mats, coeffs, eps)


def tt_displacement_fictitious(gk, sinhthetak, eps: float, pp=pp):
    r"""
    TT representation of sum_k gk[k] sinh(theta_k) (tilde a_k + tilde a_k^\dagger)
    """
    N = pp.DOF_N
    D = build_displacement_local(pp.occ)
    local_mats = [D] * N
    coeffs = gk * sinhthetak
    return tt_sum_local_operators(N, pp.occ, local_mats, coeffs, eps)


def tt_lift_physical_with_fictitious(tt_boson, left_op, eps: float, pp=pp):
    """
    Construct (left_op (x) tt_boson (x) I).
    """
    tt_left = tt_matrix(left_op)
    tt = tt_kron(tt_left, tt_boson)
    tt = tt_kron(tt, tt_eye(pp.DOF_N, pp.occ))
    return tt


def tt_lift_fictitious_with_physical(tt_boson, left_op, eps: float, pp=pp):
    """
    Construct (left_op x I x tt_boson).
    """
    tt_left = tt_matrix(left_op)
    tt = tt_kron(tt_left, tt_eye(pp.DOF_N, pp.occ))
    tt = tt_kron(tt, tt_boson)
    return tt


def tt_hamiltonian(eps: float = 1e-14, pp=pp):
    """
    Build -iH for the TFD spin-boson model using modular building blocks.

    Returns
    -------
    MPO (MPS-like TT object)
    """
    # --- parameters ---
    freq, ck, gk, thetak, sinhthetak, coshthetak = discretize_ohmic(pp.DOF_N)

    # --- electronic part ---
    He = build_electronic_hamiltonian(pp.EPSILON, pp.GAMMA_DA)
    tt_He = tt_matrix(He)
    tt_He = tt_embed_electronic(tt_He, pp.DOF_N, pp.occ)

    # --- physical and fictitious number operators ---
    tt_num_physical = tt_number_operator_physical(freq, eps)
    tt_Ie = tt_matrix(np.eye(2, dtype=np.complex128))

    tt_systemnumoc = tt_kron(tt_Ie, tt_num_physical)
    tt_systemnumoc = tt_kron(tt_systemnumoc, tt_eye(pp.DOF_N, pp.occ))

    tt_tildenumoc = tt_kron(tt_Ie, tt_eye(pp.DOF_N, pp.occ))
    tt_tildenumoc = tt_kron(tt_tildenumoc, tt_num_physical)

    # --- displacement operators ---
    tt_energy = tt_displacement_physical(gk, coshthetak, eps)
    tt_systemenergy = tt_kron(tt_matrix(np.array([[1, 0], [0, -1]],
                                                 dtype=np.complex128)),
                              tt_energy)
    tt_systemenergy = tt_kron(tt_systemenergy, tt_eye(pp.DOF_N, pp.occ))

    tt_tilenergy = tt_displacement_fictitious(gk, sinhthetak, eps)
    tt_tildeenergy = tt_kron(tt_matrix(np.array([[1, 0], [0, -1]],
                                                dtype=np.complex128)),
                             tt_eye(pp.DOF_N, pp.occ))
    tt_tildeenergy = tt_kron(tt_tildeenergy, tt_tilenergy)

    # --- assemble H ---
    H = add_tensor(tt_He, tt_systemnumoc, small=eps)
    H = add_tensor(H, tt_tildenumoc, coeff=-1.0, small=eps)
    H = add_tensor(H, tt_systemenergy, coeff=1.0, small=eps)
    H = add_tensor(H, tt_tildeenergy, coeff=1.0, small=eps)

    # fold -i into the first core
    H.nodes[0] *= -1j

    # convert to MPO and truncate
    A = MPS2MPO(H).truncation(small=eps)
    return A


def tt_ksl_propagator(
    y0: Any,
    A: Any,
    update_type: str = "rk4",
    rk4slices: int = 1,
    mmax: int = 4,
    RDO_arr_bench: np.ndarray | None = None,
    property_fn: Callable[[Any], np.ndarray] = cal_property,
    verbose: bool = True,
    show_steptime: bool = False,
    copy_state: bool = False,
    pp=pp,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Perform TT-TFD time propagation with a given initial state and Hamiltonian.

    Parameters
    ----------
    y0
        Initial TT/MPS state. If ``copy_state`` is False, this object is
        updated in-place by the propagator.
    A
        TT/MPO representing the (possibly non-Hermitian) generator, e.g. -iH.
    update_type : {"rk4", "krylov"}, optional
        Local time-stepper used in tdvp1site. Default is "rk4".
    rk4slices : int, optional
        Number of sub-slices for RK4 integration. Ignored for "krylov".
    mmax : int, optional
        Krylov subspace dimension for "krylov" updates. Default is 4.
    RDO_arr_bench : np.ndarray, optional
        Optional benchmark reduced density operator array of shape
        (TIME_STEPS, DOF_E_SQ). If provided, each step's RDO is compared
        with this reference via compare_diff.
    property_fn : callable, optional
        Function mapping the TT/MPS state to a (flattened) RDO array of shape
        (DOF_E_SQ,). Default is cal_property.
    verbose : bool, optional
        If True, print high-level progress information.
    show_steptime : bool, optional
        If True, print wall-clock time for each TDVP step.
    copy_state : bool, optional
        If True, work on a copy of `y0` instead of modifying it in-place.

    Returns
    -------
    t : np.ndarray
        1D array of simulation times of length pp.TIME_STEPS.
    RDO_arr : np.ndarray
        2D array of reduced density matrices over time with shape
        (pp.TIME_STEPS, pp.DOF_E_SQ).
    """
    n_steps = pp.TIME_STEPS
    dt = pp.DT

    # Optional copy so caller can keep original y0
    if copy_state and hasattr(y0, "copy"):
        y = y0.copy()
    else:
        y = y0

    RDO_arr = np.zeros((n_steps, pp.DOF_E_SQ), dtype=np.complex128)
    t = np.linspace(0.0, (n_steps - 1) * dt, n_steps, dtype=float)

    start_time = time.time()
    if verbose:
        print("Start propagation")
        print(f"  steps = {n_steps}, dt = {dt}, update_type = {update_type}")

    for ii, ti in enumerate(t):
        if verbose:
            print(f"Step {ii:6d}, t = {ti:.6f}")

        step_t0 = time.time()

        # TDVP one-site update
        y = tdvp1site(
            y,
            A,
            dt,
            update_type=update_type,
            mmax=mmax,
            rk4slices=rk4slices,
        )

        # Reduced density operator (or whatever property_fn returns)
        RDO_arr[ii] = property_fn(y)

        # Optional benchmark comparison
        if RDO_arr_bench is not None:
            compare_diff(RDO_arr[ii], RDO_arr_bench[ii])

        if show_steptime:
            print("  time for tdvp:", time.time() - step_t0)

    if verbose:
        print("\tTotal propagation time:", time.time() - start_time)

    return t, RDO_arr



