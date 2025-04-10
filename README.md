[![License: GNU AGPL v3](https://img.shields.io/badge/License-GNU_AGPL_v3-lightgrey.svg)](LICENCE)
[![Static Badge](https://img.shields.io/badge/CQDMQD-00268d?style=flat&logoColor=00268d&label=NSF&labelColor=00268d&color=00268d&link=https%3A%2F%2Fcqdmqd.yale.edu%2F)](https://cqdmqd.yale.edu/)


# QFlux - A Quantum Computer Dynamics Package

This repository contains various protocols for performing quantum dynamics simulations with quantum devices. Each submodule contains object implementations for these protocols as demonstrated in a publication, as well as comprehensive tutorial notebooks designed to help users understand, implement and build upon various simulation techniques for studying quantum dynamics using quantum computer frameworks. Each tutorial is provided in Python, using Jupyter Notebooks to offer detailed explanations in both markdown and code comments.


## Table of Contents

1. [Getting Started](#start)
   - [Documentation](#docs)
   - [Notebooks For Tutorial Manuscript](#notebooks)
   - [Additional Repositories](#repos)
2. [Contribution Guidelines](#contribute)
3. [Citation](#citation)
4. [License](#license)
5. [Acknowledgements](#acknowledgement)


## Getting Started <a name="start"></a>

Simply select a notebook and execute them locally or in google collab. Necessary dependencies will be installed using `pip`.

If using uv through the commandline, use the following syntax to create and activate a virtual environment:

```bash
uv venv
source .venv/bin/activate
```

The necessary packages, including development, can be installed as follows:

```bash
uv pip install -e ".[dev]"
```

### Documentation <a name="docs"></a>

Documentation for QFlux, illustrating its features and representative examples, is available at the following page:

https://qflux.batistalab.com/

### Notebooks For Tutorial Manuscript <a name="notebooks"></a>

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/tutorial/Part_I.ipynb) | Part I - Closed Quantum Dynamics and Simulation Protocols

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/tutorial/Part_I_appendixA.ipynb) | Part I - Appendix A: Variational Methods

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/tutorial/Part_I_appendixB.ipynb) | Part I - Appendix B: Bosonic Simulation

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/tutorial/Part_II.ipynb) | Part II - Open Quantum Dynamics

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/tutorial/Part_III_B.ipynb) | Part III - Variational Quantum Trajectory Dynamics for FMO

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/tutorial/Part_III_A.ipynb) | Part III - Variational Quantum Dynamics for Amplitude Damping

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/tutorial/Part_IV.ipynb) | Part IV - Generalized Quantum Master Equation Dynamics

### Contribution Guidelines <a name="contribute"></a>

To contribute to the repository, clone the qflux repository and create a new branch:

```bash
git clone git@github.com:batistagroup/qflux.git
git pull
git checkout -b part_II_docs_spinchain
```

For example to contribute to the documentation, written in markdown (.md), edit the markdown in some form and make it complete (vim spin_chain.md for example).
Once complete, commit changes to file: 

```bash
git add spin_chain.md
git commit -m 'DOCS: Added docs on the spin chain example with Lindblad'
git push
```

Generate a new pull request through github and assign to a tentative reviewer.

### Additional Repositories <a name="repos"></a>

[![Static Badge](https://img.shields.io/badge/Open_in_Github-181717.svg?&logo=github&logoColor=white)](https://github.com/dcabral00/qc_spin_tutorial) | Spin Chain Tutorial Repository 

[![Static Badge](https://img.shields.io/badge/Open_in_Github-181717.svg?&logo=github&logoColor=white)](https://github.com/saurabhshivpuje/QMAD) | QMultiAdapt Repository

[![Static Badge](https://img.shields.io/badge/Open_in_Github-181717.svg?&logo=github&logoColor=white)](https://github.com/XiaohanDan97/CCI_PartIII_GQME) | GQME Tutorial Repository


## Citation <a name="citation"></a>

Please cite the preprint of our work when using this code until the journal version becomes available.


## Licensing <a name="license"></a>

Each notebook or repository might have its own licensing. Please refer to the individual README files and notebooks within each directory for specific licensing information.

## Acknowledgement <a name="acknowledgement"></a>

We acknowledge the financial support of the National Science Foundation under award number 2124511, CCI Phase I: NSF Center for Quantum Dynamics on Modular Quantum Devices (CQD-MQD).
