# Welcome to the QFlux Documentation

![Logo](img/qflux-logo.png)

This is a Python package containing various protocols for performing quantum dynamics simulations with quantum devices. Each submodule contains object-oriented implementations for these protocols as demonstrated in our publication, as well as comprehensive tutorial notebooks designed to help users understand, implement and build upon various simulation techniques for studying quantum dynamics using quantum computer frameworks. Each tutorial is provided in Python, using Jupyter Notebooks to offer detailed explanations in both markdown and code comments.

## Installation

To get started with `qflux`, you only need to install via `pip`:

```
pip install qflux
```

This will install all necessary dependencies and set you up to check out our examples: 

- **Closed Systems**
    - [Anatomy of a Dynamics Simulation: Advanced Use Cases](basics.md)
    - [Example: Quantum Harmonic Oscillator](qho_example.md)
    - [Example: Adenine-Thymine Base Pair](AT_basepair.md)
    - [Example: Spin Chain](spinchain.md)
    - [Example: Dynamics for an Arbitrary Hamiltonian](arbitrary_evo.md)
    - [API Documentation](cs_api.md)

- **Open Systems**
    - [Open System Dynamics Overview](basics.md)
    - [Spin Chain Demo](spinchainOpen.md)
    - [Spin 1/2 Demo](spinhalfOpen.md)
    - [Double Well Demo](DoubleWellOpen.md)
    - [API Documentation](os_api.md)

- **Variational Methods**
    - [Example: Variational Quantum Time Evolution](varQTE.md)
    - [Example: Unrestricted Adaptive Variational Quantum Dynamics in an Amplitude Damping Channel](Vectorized_Adaptive.md)
    - [Example: Stochastic Schrodinger Equation for Open Systems](trajectory_FMO.md)

- **Generalized Quantum Master Equation (GQME)**
    - [Introduction to GQME](What_is_GQME.md)
    - [Introduction to TT-TFD](What_is_TTTFD.md)
    - [GQME for Spin-Boson model](spin_boson_GQME.md)
    - [Quantum Algorithms of GQME](quantum_GQME_dilation.md)

## Questions, Issues, and Feature Requests

This package and its accompanying documentation are part of an actively maintained research software project.
If you encounter unexpected behavior, identify a reproducible bug, or have questions not addressed in the documentation, please open a new issue using the [issue template](ISSUE_TEMPLATE.md). 
We also welcome suggestions for new features or improvements that could enhance the functionality, usability, or scientific scope of the codebase.

## Development Setup

If you're interested in contributing to `qflux`, you should checkout the source code and our [Contribution Guide](CONTRIBUTING.md). If you are unsure where to start, consider browsing the list of open issues or proposed enhancements to identify areas where your expertise could be most impactful.

### Development Installation

The first step to prepare your development environment is to clone the github repo. 

Once you've cloned the github repo, you should go to that directory and begin setting up your development environment with the following commands: 

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

