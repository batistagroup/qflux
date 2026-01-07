# Welcome to the QFlux Documentation

![Logo](img/qflux-logo.png)

This is a Python package containing various protocols for performing quantum dynamics simulations with quantum devices. Each submodule contains object-oriented implementations for these protocols as demonstrated in our publication, as well as comprehensive tutorial notebooks designed to help users understand, implement and build upon various simulation techniques for studying quantum dynamics using quantum computer frameworks. Each tutorial is provided as a Jupyter Notebook in Python to offer detailed explanations in both markdown and code comments.

## Installation

To get started with `qflux`, we recommend that you set-up a virtual environment with Python 3.10+. For example, you could create a conda environment: 

```bash
conda create -n qflux_env python=3.12
conda activate qflux_env
```

Once this is completed, qflux can be installed with `pip`:

```
pip install qflux
```

This will install all necessary dependencies and you are now ready to get started!

### Getting Started

You can learn how to use `qflux` through the tutorial-format examples on this documentation website:

- **Closed Systems**
    - [Anatomy of a Dynamics Simulation: Advanced Use Cases](qflux/Closed_Systems/basics.md)
    - [Example: Quantum Harmonic Oscillator](qflux/Closed_Systems/qho_example.md)
    - [Example: Adenine-Thymine Base Pair](qflux/Closed_Systems/AT_basepair.md)
    - [Example: Spin Chain](qflux/Closed_Systems/spinchain.md)
    - [Example: Dynamics for an Arbitrary Hamiltonian](qflux/Closed_Systems/arbitrary_evo.md)
    - [API Documentation](qflux/Closed_Systems/cs_api.md)

- **Open Systems**
    - [Open System Dynamics Overview](qflux/Open_Systems/basics.md)
    - [Spin Chain Demo](qflux/Open_Systems/spinchainOpen.md)
    - [Spin 1/2 Demo](qflux/Open_Systems/spinhalfOpen.md)
    - [Double Well Demo](qflux/Open_Systems/DoubleWellOpen.md)
    - [API Documentation](qflux/Open_Systems/os_api.md)

- **Variational Methods**
    - [Example: Variational Quantum Time Evolution](qflux/Variational_Methods/varQTE.md)
    - [Example: Unrestricted Adaptive Variational Quantum Dynamics in an Amplitude Damping Channel](qflux/Variational_Methods/Vectorized_Adaptive.md)
    - [Example: Stochastic Schrodinger Equation for Open Systems](qflux/Variational_Methods/trajectory_FMO.md)

- **Generalized Quantum Master Equation (GQME)**
    - [Introduction to GQME](qflux/GQME/What_is_GQME.md)
    - [Introduction to TT-TFD](qflux/GQME/What_is_TTTFD.md)
    - [GQME for Spin-Boson model](qflux/GQME/spin_boson_GQME.md)
    - [Quantum Algorithms of GQME](qflux/GQME/quantum_GQME_dilation.md)


You can also open our interactive tutorial notebooks in Google Colab: 

- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg) Closed Systems Demonstration](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/manuscript/Extended-Version-Notebooks/Part_I_qflux_snippets.ipynb)

- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg) Open Systems Demonstration](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/manuscript/Extended-Version-Notebooks/Part_II_qflux_snippets.ipynb)

- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg) Variational Methods Demonstration](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/manuscript/Extended-Version-Notebooks/Part_III_qflux_snippets.ipynb)

- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg) GQME Demonstration](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/manuscript/Extended-Version-Notebooks/Part_IV_qflux_snippets.ipynb)


We also have companion notebooks for the pre-print manuscripts available:

- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg) Classical Foundations for Quantum Dynamics Simulations - Building Intuition and Computational Workflows](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/manuscript/Part_I_JCE.ipynb)

- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg) Quantum Circuit Implementations of Molecular Dynamics - Closed Quantum Systems](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/manuscript/Part_II_JCE.ipynb)

- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg) Quantum Circuit Implementations of Molecular Dynamics - State Initialization and Unitary Decomposition](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/manuscript/Part_III_JCE.ipynb)

- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg) Dilation Method for Open Quantum Systems](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/manuscript/Part_IV_JCE.ipynb)

- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg) QMAD - A Module for Adaptive Variational Quantum Algorithms](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/manuscript/Part_V_JCE.ipynb)

- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg) The Generalized Quantum Master Equation](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/manuscript/Part_VI_JCE.ipynb)


## Questions, Issues, and Feature Requests

This package and its accompanying documentation are part of an actively maintained research software project.
If you encounter unexpected behavior, identify a reproducible bug, or have questions not addressed in the documentation, please open a new issue using the [issue template](ISSUE_TEMPLATE.md). 
We also welcome suggestions for new features or improvements that could enhance the functionality, usability, or scientific scope of the codebase.

## Development Setup

If you're interested in contributing to `qflux`, you should checkout the source code and our [Contribution Guide](CONTRIBUTING.md). If you are unsure where to start, consider browsing the list of open issues or proposed enhancements to identify areas where your expertise could be most impactful.

### Development Installation

The first step to prepare your development environment is to clone the [GitHub repo](https://github.com/batistagroup/qflux) using your preferred method (web/command line interface). For additional information on how to clone a GitHub repo, you can checkout the [GitHub Docs Page](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository). 

Once you've cloned the GitHub repo, you should go to that directory and begin setting up your development environment with the following commands: 

```bash
cd qflux-master/
```

This project uses `uv` for fast and reliable Python package management. To set up your development environment:

```bash
# Create and activate a virtual environment
uv venv
source .venv/bin/activate

# Install the package and all development dependencies
uv pip install -e ".[dev]"
# Initiate pre-commit checks
pre-commit install
uv sync
```

This will install all necessary dependencies, including development tools like pre-commit hooks, testing frameworks, and documentation generators.

### Documentation

This project uses MkDocs with the Material theme for documentation. To work with the documentation locally:

1. Make sure you have all development dependencies installed
2. Run the documentation server:

   ```bash
   mkdocs serve
   ```

3. Open your browser and navigate to `http://127.0.0.1:8000`

The documentation will automatically reload when you make changes to the markdown files.

### Code Quality Tools

We use pre-commit hooks to ensure code quality and consistency. The following tools are configured in `.pre-commit-config.yaml`:

- **Ruff**: A fast Python linter and formatter
  - Runs linting checks with auto-fix capability
  - Handles code formatting

After installing the development dependencies (as described in the Installation section), enable the pre-commit hooks by running:

```bash
pre-commit install
```

Now the hooks will run automatically on every commit, ensuring code quality and consistency.

### Managing Dependencies

This project uses `uv` for fast and reliable dependency management. Here's how to manage your dependencies:

#### Adding New Dependencies

To add a new package dependency:

```bash
uv add package_name
# Add a development dependency
uv add --dev package_name
```

This will:

1. Install the package in your virtual environment
2. Update your `pyproject.toml` with the new dependency
3. Update the `uv.lock` file with exact versions

#### Synchronizing Dependencies

If you pull changes that include new dependencies or switch branches, synchronize your environment:

```bash
uv sync
```

This ensures your virtual environment exactly matches the dependencies specified in the lock file, removing any packages you don't need and installing any that are missing.

### Writing Documentation

This project follows a structured approach to documentation. Each module should have its own markdown file in the `docs/flux/` directory, organized by topic. Documentation files might include:

1. **Overview**: A brief description of the module's purpose and key features
2. **Concepts**: Explanation of important concepts and design decisions
3. **Examples**: Code examples showing common usage patterns
4. **Source Code**: Auto-generated documentation from source code annotations

## Citing QFlux

Please cite the preprint of our work when using this code until the journal version becomes available. (We will add pre-formatted citation here once the pre-print goes live.)

## Acknowledgement of Funding

We acknowledge the financial support of the National Science Foundation under award number 2124511, CCI Phase I: NSF Center for Quantum Dynamics on Modular Quantum Devices (CQD-MQD).
