# What is GQME?
In realistic physical and chemical processes, the system of interest is typically in contact with a thermal environment. This environment involves an enormous number of degrees of freedom—for example, think of the countless solvent molecules in a solution, on the order of Avogadro’s number—making it infeasible to solve the full quantum dynamics directly. However, we are not interested in the environment for its own sake, but only insofar as it influences the dynamics of the system.

One natural strategy is to derive an equation of motion for the system alone, by effectively tracing out the environmental degrees of freedom. Such closed-form equations of motion for the reduced system are known as quantum master equations.

The **Generalized Quantum Master Equation (GQME)** is a formally exact framework for simulating reduced quantum dynamics by projecting the full Liouville-von Neumann equation onto a relevant subsystem.

In this documentation, we focus on the well-known Nakajima-Zwanzig GQME ([Nakajima, S. Prog. Theor. Phys. 1958](https://academic.oup.com/ptp/article/20/6/948/1930693), [Zwanzig, R. J.
Chem. Phys. 1960](https://doi.org/10.1063/1.1731409)):
$$
\frac{d}{dt}\hat{\sigma}(t) = -\frac{i}{\hbar}\langle \mathcal{L}\rangle_n^0\hat{\sigma}(t) - \int_0^t d\tau\, \mathcal{K}(\tau)\hat{\sigma}(t - \tau) + \mathcal{I}(t)
$$

- $\hat{\sigma}$: reduced density matrix of the subsystem of interest.
- $\langle \mathcal{L}\rangle_n^0$: the projected Liouvillian.
- $\mathcal{K}(t)$: **memory kernel**, contains all the memory effects from the bath.
- $\mathcal{I}(t)$: inhomogeneous term (often vanishes for factorized initial states).

In the following, we briefly outline its derivation and the meaning of each term. For a detailed derivation, see [this paper](https://doi.org/10.1021/acs.jctc.2c00892).

## Derivation of the Nakajima–Zwanzig Equation

We begin with the quantum Liouville–von Neumann equation, which governs the time evolution of the total system density operator:

$$
\frac{d}{dt} \hat{\rho}(t) = -\frac{i}{\hbar} [\hat{H}, \hat{\rho}(t)] \equiv -i \mathcal{L} \hat{\rho}(t),
$$

where $\mathcal{L} \equiv \hbar^{-1}[\hat{H}, \cdot]$ is the Liouvillian superoperator.

Now suppose we are only interested in a particular subsystem $S$ while treating the rest of the system as the surrounding environment $B$. The total Hamiltonian can then be partitioned as:
$$
\hat{H} = \hat{H}_S + \hat{H}_B + \hat{H}_I
$$
where $\hat{H}_S$ and $\hat{H}_B$ are the system and bath Hamiltonians, respectively, and $\hat{H}_I$ denotes their interaction.


To derive an equation of motion for the subsystem alone, we introduce a projection operator $\mathcal{P}$  that projects any operator onto the relevant subspace (i.e., the part we care about). Its complement is given by $\mathcal{Q} = 1 - \mathcal{P}$.

Applying these projectors to both sides of the Liouville equation yields a pair of coupled equations:

$$
\frac{d}{dt} \mathcal{P} \hat{\rho}(t) = -i \mathcal{P} \mathcal{L} \mathcal{P} \hat{\rho}(t) - i \mathcal{P} \mathcal{L} \mathcal{Q} \hat{\rho}(t),
$$

$$
\frac{d}{dt} \mathcal{Q} \hat{\rho}(t) = -i \mathcal{Q} \mathcal{L} \mathcal{P} \hat{\rho}(t) - i \mathcal{Q} \mathcal{L} \mathcal{Q} \hat{\rho}(t).
$$

These equations describe how the projected components $\mathcal{P} \hat{\rho}(t)$ and $\mathcal{Q} \hat{\rho}(t)$ evolve and interact over time.

Formally solving the second equation gives:
$$
\mathcal{Q} \hat{\rho}(t) = e^{-i \mathcal{Q} \mathcal{L} (t - t_0)} \mathcal{Q} \hat{\rho}(t_0)
- i \int_{t_0}^{t} d\tau\, e^{-i \mathcal{Q} \mathcal{L} (t - \tau)} \mathcal{Q} \mathcal{L} \mathcal{P} \hat{\rho}(\tau)
$$

Substituting this solution back into the first equation yields a formally closed equation for the relevant component $\mathcal{P} \hat{\rho}(t)$:

$$
\frac{d}{dt} \mathcal{P} \hat{\rho}(t) = -i \mathcal{P} \mathcal{L}\mathcal{P} \hat{\rho}(t)
- \int_{t_0}^{t} d\tau\, \mathcal{P} \mathcal{L} e^{-i \mathcal{Q} \mathcal{L}(t - \tau)} \mathcal{Q} \mathcal{L} \mathcal{P} \hat{\rho}(\tau)
- i \mathcal{P} \mathcal{L} e^{-i \mathcal{Q} \mathcal{L}(t - t_0)} \mathcal{Q} \hat{\rho}(t_0)
$$

This is the **Nakajima-Zwanzig GQME** with:
- $\hat{\sigma}(t) \equiv \mathcal{P} \hat{\rho}(t)$  the reduced density matrix,
- $\mathcal{K}(t)\equiv \mathcal{P} \mathcal{L} e^{-i \mathcal{Q} \mathcal{L}(t - \tau)} \mathcal{Q} \mathcal{L} \mathcal{P}$ the memory kernel,
- $\hbar^{-1}\langle \mathcal{L}\rangle_n^0\equiv \mathcal{P} \mathcal{L} \mathcal{P}$ the projected Liouvillian, which governs the system's unitary evolution in the absence of coupling to the bath,
- $\mathcal{I}(t) \equiv - i \mathcal{P} \mathcal{L} e^{-i \mathcal{Q} \mathcal{L}(t - t_0)} \mathcal{Q} \hat{\rho}(t_0)$ the inhomogeneous term.

This equation serves as the formal foundation for computing reduced system dynamics with full inclusion of environment-induced memory effects.


## GQME for molecular systems

For molecular systems with an overall Hamiltonian of the following form:
$$
\hat{H} = \sum_{j=1}^{N_e} \hat{H}_j(\hat{\mathbf{R}},\hat{\mathbf{P}}) |j\rangle\langle j|+ \sum_{\substack{j,k = 1 \\ k \neq j}}^{N_e} \hat{V}_{jk}(\hat{\mathbf{R}})|j\rangle\langle k|
$$
Here, $\hat{H}_j (\hat{\mathbf{R}},\hat{\mathbf{P}})= \hat{\mathbf{P}}^2/2+ V_j(\hat{\mathbf{R}})$ is the nuclear Hamiltonian when the system is in the electronic state $| j \rangle$, with the index $j$ running over the $N_e$ electronic states;
$\{ \hat{V}_{jk} (\hat{\mathbf{R}})| j \neq k \}$ are coupling terms between electronic states;
and $\hat{\mathbf{R}} = \{\hat{R}_1,\hat{R}_2,...,\hat{R}_{N_n} \}$ and $\hat{\mathbf{P}} = \{\hat{P}_1,\hat{P}_2,...,\hat{P}_{N_n} \}$ are the mass-weighted position and momentum operators of the $N_n$ nuclear degrees of freedom (DOFs).

We choose the total initial state to be a product state:
$$
\hat{\rho} (0) = \hat{\rho}_n (0) \otimes \hat{\sigma} (0)
$$
Here, $\hat{\rho}_n (0)$ is the initial density operator for the nuclear (bath) degrees of freedom, which is always taken to be in its own thermal equilibrium. Similarly, $\hat{\sigma} (0) = \text{Tr}_n \{ \hat{\rho} (0)\}$ is the reduced density operator that describes the initial state of the electronic DOFs. Here, ${\rm Tr}_n\{\cdot \}$ denotes the partial trace over nuclear (bath) degrees of freedom.

Now define the projection operator as:
$$\mathcal{P} \hat{\rho}(t) = \hat{\rho}_n(0) \otimes \text{Tr}_n\{\hat{\rho}(t)\}$$

This choice ensures that $\mathcal{P} \hat{\rho}(0) = \hat{\rho}(0)$, which means the inhomogeneous term $\mathcal{I}(t)$ in the GQME vanishes.

With this projection operator, the definitions of the projected Liouvillian and the memory kernel become:

$$
\langle \mathcal{L}\rangle_n^0 = \text{Tr}_n \left\{ \hat{\rho}_n (0) {\cal L} \right\}
$$

$$
\mathcal{K}(t) = \frac{1}{\hbar^2}\text{Tr}_n \Big\{ \mathcal{L} e^{-i \mathcal{QL} \tau / \hbar}\mathcal{QL} \hat{\rho}_n (0) \Big\}
$$



It is important to note that all effects of the environment are fully encoded in the memory kernel. In other words, the challenge of simulating an open quantum system is reduced to accurately computing the memory kernel. Once the memory kernel is known, the GQME can be directly solved to obtain the time evolution of the reduced density matrix.


### Numerical calculation of the memory kernel
We now introduce a numerical approach for evaluating the memory kernel. Using the identity for the $e^{-i \mathcal{QL} \tau / \hbar}$ [(Qiang Shi 2003)](https://doi.org/10.1063/1.1624830),[(Ming-Liang Zhang 2006)](https://doi.org/10.1063/1.2218342),

$$
e^{-i \mathcal{Q}\mathcal{L}t/\hbar}
= e^{-i \mathcal{L}t/\hbar}
+ \frac{i}{\hbar}\int_{0}^{t} \mathrm{d}\tau \,
e^{-i \mathcal{L}(t - \tau)/\hbar} \mathcal{P}\mathcal{L} e^{-i \mathcal{Q}\mathcal{L}\tau/\hbar}
$$

The memory kernel becomes
$$
\mathcal{K}(t)
= i\dot{\mathcal{F}}(t)
- \frac{1}{\hbar}\mathcal{F}(t)\langle \mathcal{L}\rangle_n^0
+ i\int_{0}^{t} d\tau \, \mathcal{F}(t - \tau)\mathcal{K}(\tau)
$$
with the projection-free inputs (PFIs) $\mathcal{F}(t)$ and $\dot{\mathcal{F}}(t) $ defined as
$$
\mathcal{F}(t)
= \frac{1}{\hbar}\mathrm{Tr}_n\bigl[\mathcal{L} e^{-i\mathcal{L}t/\hbar}\hat{\rho}_n(0)\bigr]
$$
$$
\dot{\mathcal{F}}(t)
= -\frac{i}{\hbar^2}\mathrm{Tr}_n\bigl[\mathcal{L} e^{-i\mathcal{L}t/\hbar}\mathcal{L}\hat{\rho}_n(0)\bigr].
$$

The projection-free inputs can be related to the dynamics of the reduced density operator.
Note that the reduce density operator $\hat{\sigma}(t)$ can be formally written as
$$
\hat{\sigma}(t)  = \mathcal{G}(t) \hat{\sigma}(0) = \mathrm{Tr}_n [e^{-i\mathcal{L}t/\hbar}\hat{\rho}_n(0)] \hat{\sigma}(0)
$$
therefore, $\mathcal{F}(t) = i\dot{\mathcal{G}}(t)$. $\mathcal{F}(t)$ and $\dot{\mathcal{F}}(t) $ can be obtained through taking time-derivatives of $\mathcal{G}(t)$.

The propagator $\mathcal{G}(t)$, is a super-operator with the matrix
element $\mathcal{G}_{jk,lm}(t)$, which can be defined by starting from initial state $|l⟩⟨m| ⊗ \hat{\rho}_n(0)$, measure the $\sigma_{jk}(t)$ at time $t$.

These quantities can be computed using various existing numerical methods, such as the numerically exact tensor-train thermo-field dynamics (TT-TFD) approach introduced in  [TT-TFD](What_is_TTTFD.md). We will demonstrate the practical solution of the GQME using `qflux` in [GQME for Spin-Boson model](spin_boson_GQME.md).

# Summary
This example introduces the Generalized Quantum Master Equation (GQME) framework, which is an important approach for simulating open quantum system dynamics:

* The Nakajima–Zwanzig formalism is introduced as the theoretical foundation for deriving the GQME.

* For chemically relevant molecular systems, we provide explicit definitions of each term in the GQME.

* The memory kernel expression is expanded, illustrating how it can be computed numerically using various methods.

The GQME framework thus provides a rigorous and versatile foundation for investigating quantum dynamics in systems coupled to realistic environments.
