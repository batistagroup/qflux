# Library containing RIDA, TREX, and ZNE Mitigation functions

import numpy as np
import math
import matplotlib.pyplot as plt

from qiskit import QuantumCircuit, transpile
from qiskit.circuit import ClassicalRegister
from qiskit.quantum_info import Statevector
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error, ReadoutError

#-------------------------------------------------------
# Utilities
# create a random ansatz
def random_ansatz(n_qubits=5, depth=6, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    qc = QuantumCircuit(n_qubits, 1)
    for _ in range(depth):
        for q in range(n_qubits):
            qc.ry(2*np.pi*rng.random(), q)
            qc.rz(2*np.pi*rng.random(), q)
        for q in range(n_qubits):
            qc.cx(q, (q+1) % n_qubits)
    qc.measure(0, 0)
    return qc

# true expectation value of <Z> on q0
def true_expectation_z(qc):
    n = qc.num_qubits
    base = QuantumCircuit(n)
    for instr in qc.data:      # strip measurements
        inst, qargs, cargs = instr.operation, instr.qubits, instr.clbits
        if inst.name != 'measure':
            base.append(inst, qargs, cargs)
    psi = Statevector.from_instruction(base)
    probs = psi.probabilities()
    z = 0.0
    for idx, p in enumerate(probs):
        bit0 = (idx & 1)# 1 is ...0001 in binary so this extracts the Least Signifcant Bit(q0) for each bitstring
        z += (1.0 if bit0==0 else -1.0)*p
    return z

def run_counts(qc, shots, backend, noise_model):
    sim = AerSimulator(noise_model=noise_model)
    tqc = transpile(qc, backend=sim, optimization_level=0)
    return sim.run(tqc, shots=shots).result().get_counts()

#--------------------------------------------------------
# A. Zero-Noise Extrapolation (ZNE)

# Noise model with depolarizing and readout noise
def make_noise_model(p1=2.5e-4, p2=2.0e-3, r_meas=7.5e-3, scale=1.0):
    """Depolarizing + readout noise. ZNE scaling via p' = 1 - (1-p)**scale."""

    # Repeat the noisy process "scale" times
    p1s = 1.0 - (1.0 - p1)**scale
    p2s = 1.0 - (1.0 - p2)**scale
    rms = 1.0 - (1.0 - r_meas)**scale

    nm = NoiseModel()
    e1 = depolarizing_error(p1s, 1)
    for g in ['rx', 'rz', 'x', 'sx', 'ry']:
        nm.add_all_qubit_quantum_error(e1, g)

    e2 = depolarizing_error(p2s, 2)
    nm.add_all_qubit_quantum_error(e2, 'cx')

    ro = ReadoutError([[1-rms, rms],
                       [rms, 1-rms]])
    nm.add_all_qubit_readout_error(ro)
    return nm

# measure z expectation for a given qubit. Qiskit convention: q0 is the rightmost bit
def measure_z_expectation(counts, bit_index=0, shots=None):
    if shots is None:
        shots = sum(counts.values())
    z = 0.0
    for bitstring, c in counts.items():
        bit = bitstring[-1 - bit_index]  # indexing from the right
        z += (1.0 if bit == '0' else -1.0) * c
    return z / shots

def zne_exponential_fit(x1, x3, x5):
    # Robust closed-form + fallbacks
    if (x3 <= x1 <= x5) or (x5 <= x1 <= x3):
        return 0.5 * (x1 + x3)
    if (x1 <= x5 <= x3) or (x3 <= x5 <= x1):
        return (3 * x1 - x3) / 2
    if np.isclose(x1, x5) and not np.isclose(x1, x3):
        return x1
    denom = (x1 - x3)
    if np.isclose(denom, 0.0):
        return 0.5 * (x1 + x3)
    u = (x3 - x5) / denom
    if u <= 0 or not np.isfinite(u):
        return (3 * x1 - x3) / 2
    return x1 + (x1 - x3) / (u + math.sqrt(u))

# New: to expose the whole workflow 
def zne_estimate(qc, shots, backend, p1, p2, r_meas):
    """
    Estimate <Z> using Zero Noise Extrapolation.
    """

    s1 = shots // 3
    s3 = shots // 3
    s5 = shots - s1 - s3

    nm1 = make_noise_model(p1, p2, r_meas, scale=1)
    nm3 = make_noise_model(p1, p2, r_meas, scale=3)
    nm5 = make_noise_model(p1, p2, r_meas, scale=5)

    O1 = measure_z_expectation(run_counts(qc, s1, backend, nm1), 0, s1)
    O3 = measure_z_expectation(run_counts(qc, s3, backend, nm3), 0, s3)
    O5 = measure_z_expectation(run_counts(qc, s5, backend, nm5), 0, s5)

    return zne_exponential_fit(O1, O3, O5)

#--------------------------------------------------------
# B. RIDA

# Convert qubit references from template circuit to equivalent qubits in a new circuit, allowing
# gates to be copied between circuits without changing their targets
def _map_qargs_to_new_circuit(qargs, template_circ, new_circ):
    mapped = []
    for q in qargs:
        pos = template_circ.find_bit(q).index
        mapped.append(new_circ.qubits[pos])
    return mapped

# Construct a RIDA estimation circuit by randomly selecting roughly half
# of the one and two-qubit gates
def build_rida_estimation_circuit(qc, backend, rng=None):
    if rng is None:
        rng = np.random.default_rng()

    no_meas = qc.remove_final_measurements(inplace=False)
    tq = transpile(no_meas, backend=AerSimulator(), optimization_level=1)
    n = tq.num_qubits

    oneq_idxs, twoq_idxs, ops = [], [], []
    for idx, instr in enumerate(tq.data):
        inst, qargs, _ = instr.operation, instr.qubits, instr.clbits
        if inst.name in ['rx','ry','rz','x','sx']:
            oneq_idxs.append(idx)
        elif inst.name=='cx':
            twoq_idxs.append(idx)
        ops.append((inst,qargs))

    k1, k2 = len(oneq_idxs)//2, len(twoq_idxs)//2
    sel1 = set(rng.choice(oneq_idxs, k1, replace=False)) if k1>0 else set()
    sel2 = set(rng.choice(twoq_idxs, k2, replace=False)) if k2>0 else set()
    selected = [i for i in range(len(ops)) if i in sel1 or i in sel2]

    est = QuantumCircuit(n,1)
    for i in selected:
        inst,qargs = ops[i]
        est.append(inst, _map_qargs_to_new_circuit(qargs,tq,est))
    for i in reversed(selected):
        inst,qargs = ops[i]
        try: inv = inst.inverse()
        except: inv = inst
        est.append(inv, _map_qargs_to_new_circuit(qargs,tq,est))
    est.measure(0,0)
    return est, tq

# New -- for consistency with the ZNE and TREX
def rida_estimate(qc, shots, backend, noise_model, rng=None):
    """
    Estimate <Z> using Random Identity Decomposition Approximation (RIDA).
    """

    if rng is None:
        rng = np.random.default_rng()

    s_target = shots // 2
    s_est = shots - s_target

    # Run target circuit
    counts_target = run_counts(qc, s_target, backend, noise_model)
    O_target = measure_z_expectation(counts_target, 0, s_target)

    # Build and run estimation circuit
    est_circuit, _ = build_rida_estimation_circuit(qc, backend, rng=rng)
    counts_est = run_counts(est_circuit, s_est, backend, noise_model)
    O_est = measure_z_expectation(counts_est, 0, s_est)

    if np.isclose(O_est, 0.0):
        return O_target

    return O_target / O_est

#------------------------------------------------------------------
# C. TREX

# Estimate <Z> using TREX by correcting readout errors through measurement
# twirling and calibration-based normalization.
def trex_estimate(qc, shots, backend, noise_model):

    # make two circuits, in one of them apply X before measurement
    qcA = qc.copy()
    qcB = qc.remove_final_measurements(inplace=False)
    if qcB.num_clbits==0: qcB.add_bits([ClassicalRegister(1)[0]])
    qcB.x(0); qcB.measure(0,0)
    sA,sB = shots//2, shots - shots//2
    zA = measure_z_expectation(run_counts(qcA,sA,backend,noise_model),0,sA)
    zB = -measure_z_expectation(run_counts(qcB,sB,backend,noise_model),0,sB) # -1 converts it back after flipping with X gate

    # weighted average of the two estimates - to surpress assymetric readout errors
    fD1 = (sA*zA+sB*zB)/shots

    # make two circuits again with known outputs (1/-1) to estimate how strongly the measurement process attenuates <Z>
    n = qc.num_qubits
    cal  = QuantumCircuit(n,1); cal.measure(0,0)
    calB = cal.remove_final_measurements(inplace=False)
    if calB.num_clbits==0: calB.add_bits([ClassicalRegister(1)[0]])
    calB.x(0); calB.measure(0,0)
    zA = measure_z_expectation(run_counts(cal ,sA,backend,noise_model),0,sA)
    zB =-measure_z_expectation(run_counts(calB,sB,backend,noise_model),0,sB)
    fD0 = (sA*zA+sB*zB)/shots
    return fD1 if np.isclose(fD0,0) else fD1/fD0


#------------------------------------------------------------------
# D. ZNE + TREX

def zne_trex_estimate(qc, shots, backend, p1, p2, r_meas):
    s1 = shots//3; s3 = shots//3; s5 = shots - s1 - s3
    nm1 = make_noise_model(p1,p2,r_meas,scale=1)
    nm3 = make_noise_model(p1,p2,r_meas,scale=3)
    nm5 = make_noise_model(p1,p2,r_meas,scale=5)
    t1 = trex_estimate(qc,s1,backend,nm1)
    t3 = trex_estimate(qc,s3,backend,nm3)
    t5 = trex_estimate(qc,s5,backend,nm5)
    return zne_exponential_fit(t1,t3,t5)

#----------------------------------------------------------------------------
# E. Comparison across techniques
def single_trial(shots, backend, p1, p2, r_meas, rng):

    # make a random circuit and estimate ideal <Z> on q0
    target = random_ansatz(n_qubits=5, depth=6, rng=rng)
    O_true = true_expectation_z(target)

    # make noisy circuit and estimate <Z> on q0
    nm = make_noise_model(p1, p2, r_meas, scale=1.0)
    counts = run_counts(target, shots, backend, nm)
    O_unmit = measure_z_expectation(counts, 0, shots)

    # apply the mitigation techniques
    O_rida = rida_estimate(target, shots, backend, nm, rng=rng)
    O_zne = zne_estimate(target, shots, backend, p1, p2, r_meas)
    O_zne_trex = zne_trex_estimate(target, shots, backend, p1, p2, r_meas)

    return (O_unmit - O_true)**2, (O_rida - O_true)**2, (O_zne - O_true)**2, (O_zne_trex - O_true)**2

def experiment(shots_list, num_trials=50, seed=7):
    rng = np.random.default_rng(seed)
    backend = AerSimulator()

    p1 = 2.5e-4    # 1q depol
    p2 = 2.0e-3    # 2q depol
    r_meas = 7.5e-3

    rmse_u, rmse_r, rmse_z, rmse_zt = [], [], [], []
    for s in shots_list:
        errs_u, errs_r, errs_z, errs_zt = [], [], [], []
        for _ in range(num_trials):
            eu, er, ez, ezt = single_trial(s, backend, p1, p2, r_meas, rng)
            errs_u.append(eu); errs_r.append(er); errs_z.append(ez); errs_zt.append(ezt)
        rmse_u.append(np.sqrt(np.mean(errs_u)))
        rmse_r.append(np.sqrt(np.mean(errs_r)))
        rmse_z.append(np.sqrt(np.mean(errs_z)))
        rmse_zt.append(np.sqrt(np.mean(errs_zt)))
    return np.array(rmse_u), np.array(rmse_r), np.array(rmse_z), np.array(rmse_zt)

