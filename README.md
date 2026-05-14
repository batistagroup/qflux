[![License: GNU AGPL v3](https://img.shields.io/badge/License-GNU_AGPL_v3-lightgrey.svg)](LICENSE)
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

`qflux` can be installed via `pip`: 

```bash
pip install qflux
```

To get started, one can simply select a notebook and execute them locally or in google collab. Necessary dependencies will be installed using `pip`.

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

- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg) Classical Foundations for Quantum Dynamics Simulations - Building Intuition and Computational Workflows](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/manuscript/SI_I.ipynb)

- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg) Quantum Circuit Implementations of Molecular Dynamics - Closed Quantum Systems](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/manuscript/SI_II.ipynb)

- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg) Quantum Circuit Implementations of Molecular Dynamics - State Initialization and Unitary Decomposition](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/manuscript/SI_III.ipynb)

- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg) Dilation Method for Open Quantum Systems](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/manuscript/SI_IV.ipynb)

- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg) Adaptive Variational Quantum Algorithms](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/manuscript/SI_V.ipynb)

- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg) The Generalized Quantum Master Equation](https://colab.research.google.com/github/batistagroup/qflux/blob/master/demos/manuscript/SI_VI.ipynb)

### Contribution Guidelines <a name="contribute"></a>

To contribute to the repository, follow the procedure outlined in the [Contribution Guidelines](https://github.com/batistagroup/qflux/blob/master/CONTRIBUTING.md) markdown file. 

### Additional Repositories <a name="repos"></a>

This section includes additional repositories with functionality that has been integrated within QFlux.

[![Static Badge](https://img.shields.io/badge/Open_in_Github-181717.svg?&logo=github&logoColor=white)](https://github.com/dcabral00/qc_spin_tutorial) | Spin Chain Tutorial Repository 

[![Static Badge](https://img.shields.io/badge/Open_in_Github-181717.svg?&logo=github&logoColor=white)](https://github.com/saurabhshivpuje/QMAD) | QMAD Repository

[![Static Badge](https://img.shields.io/badge/Open_in_Github-181717.svg?&logo=github&logoColor=white)](https://github.com/XiaohanDan97/CCI_PartIII_GQME) | GQME Tutorial Repository


## Citation <a name="citation"></a>

Please cite the ChemRxiv preprints of our work when using this code until the journal versions become available.

```bibtex
@article{qflux.2026_1, 
  year    = {2026}, 
  title   = {{QFlux}: Classical Foundations for Quantum Dynamics Simulation. Part I - Building Intuition and Computational Workflows}, 
  author  = {Allen, Brandon C and Dan, Xiaohan and Cabral, Delmar G A and Vu, Nam P and Cianci, Cameron and Soudackov, Alexander V and Dutta, Rishab and Kais, Sabre and Geva, Eitan and Batista, Victor S}, 
  journal = {{ChemRxiv}}, 
  doi     = {10.26434/chemrxiv.10001765/v1}
}

@article{qflux.2026_2, 
  year    = {2026}, 
  title   = {{QFlux}: Quantum Circuit Implementations of Molecular Dynamics. Part {II} - Closed Quantum Systems}, 
  author  = {Cabral, Delmar G A and Allen, Brandon C and Cianci, Cameron and Soudackov, Alexander V and Dan, Xiaohan and Vu, Nam P and Dutta, Rishab and Kais, Sabre and Geva, Eitan and Batista, Victor S}, 
  journal = {{ChemRxiv}}, 
  doi     = {10.26434/chemrxiv.10001766/v1}
}

@article{qflux.2026_3, 
  year    = {2026}, 
  title   = {{QFlux}: Quantum Circuit Implementations of Molecular Dynamics. Part {III} - State Initialization and Unitary Decomposition}, 
  author  = {Soudackov, Alexander V and Cabral, Delmar G A and Allen, Brandon C and Dan, Xiaohan and Vu, Nam P and Cianci, Cameron and Dutta, Rishab and Kais, Sabre and Geva, Eitan and Batista, Victor S}, 
  journal = {{ChemRxiv}}, 
  doi     = {10.26434/chemrxiv.10001767/v1}
}

@article{qflux.2026_4, 
  year    = {2026}, 
  title   = {{QFlux}: An Open-Source Toolkit for Quantum Dynamics Simulations on Quantum Computers. Part {IV} - Dilation Method for Open Quantum Systems}, 
  author  = {Dan, Xiaohan and Shivpuje, Saurabh and Wang, Yuchen and Cabral, Delmar G A and Allen, Brandon C and Khazaei, Pouya and Soudackov, Alexander V and Hu, Zixuan and Lyu, Ningyi and Geva, Eitan and Kais, Sabre and Batista, Victor S}, 
  journal = {{ChemRxiv}}, 
  doi     = {10.26434/chemrxiv.10001768/v1}
}

@article{qfluc.2026_5, 
  year    = {2026}, 
  title   = {{QFlux}: An Open-Source Toolkit for Quantum Dynamics Simulations on Quantum Computers. Part V - Adaptive Variational Quantum Algorithms for Open Quantum Systems}, 
  author  = {Shivpuje, Saurabh and Soudackov, Alexander V and Dan, Xiaohan and Wang, Yuchen and Allen, Brandon C and Cabral, Delmar G A and Hu, Zixuan and Lyu, Ningyi and Geva, Eitan and Batista, Victor S and Kais, Sabre}, 
  journal = {{ChemRxiv}}, 
  doi     = {10.26434/chemrxiv.10001769/v2}
}

@article{qflux.2026_6, 
  year    = {2026}, 
  title   = {{QFlux}: An Open-Source Toolkit for Quantum Dynamics Simulations on Quantum Computers. Part {VI} - The Generalized Quantum Master Equation}, 
  author  = {Dan, Xiaohan and Khazaei, Pouya and Allen, Brandon C and Lyu, Ningyi and Wilson, Callie and Mulvihill, Ellen and Wang, Yuchen and Shivpuje, Saurabh and Kais, Sabre and Batista, Victor S and Geva, Eitan}, 
  journal = {{ChemRxiv}}, 
  doi     = {10.26434/chemrxiv.10001770/v1}
}

```


## Licensing <a name="license"></a>

Each notebook or repository might have its own licensing. Please refer to the individual README files and notebooks within each directory for specific licensing information.


## Acknowledgement <a name="acknowledgement"></a>

We acknowledge the financial support of the National Science Foundation under award number 2124511, CCI Phase I: NSF Center for Quantum Dynamics on Modular Quantum Devices (CQD-MQD).
