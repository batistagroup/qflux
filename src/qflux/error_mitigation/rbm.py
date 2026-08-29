# Functions defining, used by and related to a Restricted Boltzmann Machine (RBM)

import jax
import jax.numpy as jnp
from flax import linen as nn
import numpy as np
from tqdm import tqdm

# RBM CLASSES/FUNCTIONS
#-------------------------------------------------
# A. General RBM function -- Flax/Jax version

# Utility 
def binary_to_spin(x):
    """
    Convert:0 -> +1 : 1 -> -1
    """
    return 1.0 - 2.0 * x 

# Complex RBM made using flax (so that we can calculate gradients with jax later)
class RBM(nn.Module):

    alpha: int = 1
    param_dtype: any = jnp.complex128
    kernel_init_std: float = 0.01
 
    # Activation
    def activation(self, x):
        return jnp.log(jnp.cosh(x))

    # Forward pass
    @nn.compact
    def __call__(self, x):

        # Convert binary to spin convention
        x = binary_to_spin(x)
        n_visible = x.shape[-1]
        n_hidden = self.alpha * n_visible 

        # Dense hidden layer
        dense = nn.Dense(
            features=n_hidden,
            use_bias=True,
            param_dtype=self.param_dtype,
            kernel_init=nn.initializers.normal(self.kernel_init_std),
            bias_init=nn.initializers.normal(self.kernel_init_std))

        hidden = dense(x)

        # RBM hidden layer contribution
        hidden_term = jnp.sum(self.activation(hidden), axis=-1)

        # Visible bias
        visible_bias = self.param(
            "visible_bias",
            nn.initializers.normal(
                self.kernel_init_std),
            (n_visible,),
            self.param_dtype,
        )

        visible_term = jnp.dot(x, visible_bias)
        
        # Final log-wavefunction
        return hidden_term + visible_term

#-------------------------------------------------
# B. Single-Spin demo functions

# Define psi(sigma)=exp(a * sigma) and compute probabilities in X, Y, Z bases.
def psi(a, sigma):
    return np.exp(a * sigma)

def model_probs(a):
    sig = np.array([+1, -1])

    # Z basis
    amp_z = np.array([psi(a, s) for s in sig])
    pz = np.abs(amp_z)**2
    pz /= pz.sum()

    # X basis
    amp_x = np.array([(psi(a,+1) + s*psi(a,-1))/np.sqrt(2) for s in sig])
    px = np.abs(amp_x)**2
    px /= px.sum()

    # Y basis
    amp_y = np.array([(psi(a,+1) + 1j*s*psi(a,-1))/np.sqrt(2) for s in sig])
    py = np.abs(amp_y)**2
    py /= py.sum()

    return {"X": px, "Y": py, "Z": pz}

def kl(p, q, eps=1e-12):
    return float(np.sum(p * np.log((p+eps)/(q+eps))))

#-------------------------------------------------
# C. Bell State demo functions

# Analytic RBM parameters
b = 1j * np.pi / 2
W = 1j * np.pi / 4

# List out possible configurations
spins = np.array([+1, -1])
configs = [(s1, s2) for s1 in spins for s2 in spins]

def bell_psi(s1, s2):
    """RBM wavefunction amplitude psi(sigma 1, sigma 2)."""
    return 2.0 * np.cosh(b + W * (s1 - s2))


def bell_model_probs():
    """Compute RBM probabilities in ZZ, XX, YY bases."""

    # Z tensor Z
    amp_zz = np.array([bell_psi(s1, s2) for (s1, s2) in configs])
    pzz = np.abs(amp_zz)**2
    pzz /= pzz.sum()

    # X tensor X
    amp_xx = []
    for s1, s2 in configs:
        amp = 0.0
        for sp1, sp2 in configs:
            amp += (1 + s1 * sp1) * (1 + s2 * sp2) * bell_psi(sp1, sp2)
        amp_xx.append(amp / 2.0)

    amp_xx = np.array(amp_xx)
    pxx = np.abs(amp_xx)**2
    pxx /= pxx.sum()

    # Y tensor Y
    amp_yy = []
    for s1, s2 in configs:
        amp = 0.0
        for sp1, sp2 in configs:
            amp += (1 + 1j * s1 * sp1) * (1 + 1j * s2 * sp2) * bell_psi(sp1, sp2)
        amp_yy.append(amp / 2.0)

    amp_yy = np.array(amp_yy)
    pyy = np.abs(amp_yy)**2
    pyy /= pyy.sum()

    return {"ZZ": pzz, "XX": pxx, "YY": pyy}

#-------------------------------------------------------------
# D. General RBM without Flax/Jax

class ComplexRBM:
    def __init__(self, n_visible, n_hidden, rng=None, scale=0.01):
        self.N, self.M = n_visible, n_hidden
        self.rng = np.random.default_rng() if rng is None else rng
        # complex parameters
        self.a = scale*(self.rng.standard_normal(self.N) + 1j*self.rng.standard_normal(self.N))
        self.b = scale*(self.rng.standard_normal(self.M) + 1j*self.rng.standard_normal(self.M))
        self.W = scale*(self.rng.standard_normal((self.N,self.M)) + 1j*self.rng.standard_normal((self.N,self.M)))

    # Argument of cosh in closed form definition of RBM
    def theta(self, sigma):
        return self.b + sigma @ self.W

    # log-wavefunction
    def logpsi(self, sigma):
        th = self.theta(sigma)
        return np.sum(self.a * sigma) + np.sum(np.log(2.0*np.cosh(th)))

    # log of probabilities
    def log_prob_sigma(self, sigma):
        return 2.0*np.real(self.logpsi(sigma))
    
    # Returns the wavefunction
    def psi(self, s):
        return np.exp(self.logpsi(s))

    # Sample bitstrings from the RBM
    def metropolis_samples(self, n_samples=2000, burn_in=500, thin=10):
        sigma = 2*self.rng.integers(0,2,size=self.N) - 1
        th = self.theta(sigma)
        logp = 2.0*np.real(np.sum(self.a*sigma) + np.sum(np.log(2.0*np.cosh(th))))
        S = []
        steps = burn_in + n_samples*thin
        for t in range(steps):
            i = self.rng.integers(self.N)
            sigma_new = sigma.copy(); sigma_new[i] *= -1
            th_new = th - 2.0*self.W[i,:]*sigma[i]
            logp_new = 2.0*np.real(np.sum(self.a*sigma_new) + np.sum(np.log(2.0*np.cosh(th_new))))
            if logp_new >= logp or np.log(self.rng.random()) < (logp_new - logp):
                sigma, th, logp = sigma_new, th_new, logp_new
            if t >= burn_in and ((t - burn_in) % thin == 0):
                S.append(sigma.copy())
        return np.array(S, dtype=np.int8)
   
    
# NOISE MODELS AND NOISE-AWARE FUNCTIONS
#--------------------------------------------------------------
# E. Noisy readout channel and likelihood functions

def sample_noisy_from_clean(clean_bits01, M_list, rng):
    n, N = clean_bits01.shape
    y = clean_bits01.copy()
    for i in range(N):
        p0, p1 = M_list[i][0,1], M_list[i][1,0]  # flip probabilities
        flips = rng.random(n) < np.where(clean_bits01[:, i]==0, p1, p0)
        y[:, i] ^= flips.astype(np.int64)
    return y

def log_p_y_given_sigma(y_bits01, sigma_pm1, M_list):
    # log Pr(y | sigma) under independent readout model
    x = (1 - sigma_pm1) // 2  # map {+1,-1} -> {0,1}
    N = x.shape[-1]
    logp = 0.0
    for i in range(N):
        logp += np.log(M_list[i][x[i], y_bits01[i]] + 1e-30)
    return logp

def noise_aware_loss(rbm, Y_bits01, M_list, n_mc=512):
    S = rbm.metropolis_samples(n_samples=n_mc, burn_in=300, thin=5)
    logp_sigma = np.array([rbm.log_prob_sigma(s) for s in S])
    w = np.exp(logp_sigma - np.max(logp_sigma))
    w = w / (np.sum(w) + 1e-30)
    loss = 0.0
    for y in Y_bits01:
        log_terms = np.array([log_p_y_given_sigma(y, s, M_list) + np.log(w[k]+1e-300)
                              for k,s in enumerate(S)])
        m = np.max(log_terms)
        loss += m + np.log(np.sum(np.exp(log_terms - m)) + 1e-300)
    return loss / len(Y_bits01)

#------------------------------------------------------------
# F. Noise-aware training
def train_noise_aware_rbm(Y_bits01, n_visible, n_hidden, M_list, epochs=200, lr=0.03, seed=7):
    rng = np.random.default_rng(seed)
    rbm = ComplexRBM(n_visible, n_hidden, rng=rng, scale=0.02)
    base_loss = noise_aware_loss(rbm, Y_bits01, M_list, n_mc=512)
    for ep in tqdm(range(1, epochs+1), desc="Noise aware training progress", colour='green'):
        for param in ['a','b','W']:
            T = getattr(rbm, param)
            noise = (rng.standard_normal(T.shape) + 1j*rng.standard_normal(T.shape))*0.001
            setattr(rbm, param, T + noise)
            L_pos = noise_aware_loss(rbm, Y_bits01, M_list, n_mc=256)
            gain = np.real(L_pos - base_loss)
            setattr(rbm, param, T + lr*gain*noise)
        base_loss = noise_aware_loss(rbm, Y_bits01, M_list, n_mc=256)
        if ep % 20 == 0:
            print(f"Epoch {ep:3d}  L_NA ≈ {base_loss:.4f}")
    return rbm


# LOCAL ESTIMATORS
#-------------------------------------------------------------------
# G. Pauli utilities on spin configs sigma in {+1,-1}^N

def apply_pauli_on_config(sigma, pauli_word):
    sigma_p = sigma.copy()
    phase = 1.0 + 0.0j
    for i, P in enumerate(pauli_word):
        if P == 'I':
            continue
        elif P == 'Z':
            phase *= sigma[i]  # eigenvalue ±1
        elif P == 'X':
            sigma_p[i] *= -1
        elif P == 'Y':
            phase *= (1j * sigma[i])
            sigma_p[i] *= -1
        else:
            raise ValueError("Bad Pauli")
    return sigma_p, phase


#  Local estimator for O = sum_k c_k P_k
def local_estimator_sample(rbm, sigma, pauli_terms):
    logpsi_sigma = rbm.logpsi(sigma)
    psi_sigma = np.exp(logpsi_sigma)
    Oloc = 0.0 + 0.0j
    for ck, Pk in pauli_terms:
        sigma_k, phase_k = apply_pauli_on_config(sigma, Pk)
        psi_sigma_k = rbm.psi(sigma_k)
        Oloc += ck * phase_k * (psi_sigma_k / psi_sigma)
    return Oloc

def estimate_observable_mc(rbm, pauli_terms, n_samples=50000):
    S = rbm.metropolis_samples(n_samples=n_samples, burn_in=2000, thin=5)
    vals = np.array([local_estimator_sample(rbm, s, pauli_terms) for s in S])
    mean = vals.mean()
    stderr = vals.std()/np.sqrt(len(vals))
    return mean, stderr

def raw_shot_estimator(noisy_bits01, pauli_terms):
    # Supports Pauli words with only I/Z for direct estimate from {0,1} bits.
    y_pm1 = 1 - 2*noisy_bits01  # {0,1} -> {+1,-1}
    vals = []
    for ck, Pk in pauli_terms:
        if any(P not in ('I','Z') for P in Pk):
            raise ValueError("Baseline supports I/Z only.")
        term = np.ones(noisy_bits01.shape[0])
        for q, P in enumerate(Pk):
            if P == 'Z':
                term *= y_pm1[:, q]
        vals.append(ck * term.mean())
    return np.sum(vals)