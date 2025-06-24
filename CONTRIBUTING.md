# Welcome to the QFlux Contributing Guide!

Thank you for taking the time to contribute to the **QFlux** package. Your help improves the software and expands its utility for the scientific community.

This guide contains an overview of how to open issues, make code or documentation changes, and contribute via pull requests (PRs).

---

## New Contributor Guide

If you're new to open-source or QFlux, welcome! You're encouraged to start by:

* Reading existing issues.
* Examining existing example notebooks under `demos/`.
* Looking at the `good first issue` or `help wanted` tags.
* Exploring the source code under `src/qflux`.
* Reviewing [documentation](https://qflux.batistalab.com/) in the `docs/` folder (built with [MkDocs](https://www.mkdocs.org/)).

---

## Getting Started

The general procedure for contributing to `QFlux` is: 

1. Clone the repository.
2. Install dependencies.
3. Create a branch.
4. Add your edits/modifications.
5. Create a pull request.

### Clone the repository

```bash
git clone git@github.com:batistagroup/qflux.git
cd qflux
```

### Requirements

QFlux uses [uv](https://github.com/astral-sh/uv) for dependency management. You will need:

* Python ≥3.10
* [uv](https://github.com/astral-sh/uv): Install with `pip install uv`
* [mkdocs](https://www.mkdocs.org/) for documentation
* [mkdocs-material](https://squidfunk.github.io/mkdocs-material/) documentation theme

All dependencies can be installed with the following steps (assuming you have already cloned the github repository and are within the main `qflux` directory): 

1. Install `uv`:

```bash
pip install uv
```

2. Create a virtual environment:

```bash
uv venv
```

3. Activate the virtual environment:

```bash
source .venv/bin/activate
```

4. Install dependencies:

```bash
uv  pip install -e ".[dev]"
```

You should now have a proper environment with all necessary dependencies! 

### Project structure

```
qflux/
├── src/           # Core package source code, split according to subdomain (e.g., src/qflux)
├── data/          # Input/output data files for examples
├── docs/          # Markdown documentation for use with MkDocs
├── tests/         # Unit and integration tests
├── pyproject.toml # Project configuration with uv and build info
```

---

## Issues

### Create an Issue

1. Go to the [Issues](https://github.com/batistagroup/qflux/issues) tab.
2. Use an appropriate template (e.g., *bug*, *feature request*, or *question*).
3. Be concise and clear.
4. Link to relevant code lines/commit if applicable.

### Solve an Issue

* If you'd like to address an open issue, comment on it to indicate you're working on it.
* Fork the repository and follow the development workflow below.

---

## Make Changes

### New Feature

1. Create a feature branch:

   ```bash
   git checkout -b feat/qumode_dynamics
   ```
2. Implement your code under a relevant submodule of `src/qflux/`.
3. Add or update relevant tests in `tests/`.
4. Ensure the test suite passes:

   ```bash
   uv pip install -r requirements-dev.txt
   pytest
   ```

### Documentation Updates

In this section, we will provide a general workflow for updating the documentation. We will first emphasize an important consideration prior to editing the docs. You should navigate to the directory for the module that you want to edit the documentation for. In order to determine the "routing" for these markdown files, you should consult the `mkdocs` configuration file in the top-level `qflux/` directory, named `mkdocs.yml`. In this file, you'll see a `nav:` section which  specifies the pages/routing for the documentation website. Individual pages are specified as key-value pairs where the key specifies the page title and the  value specifies the corresponding markdown file that the page is constructed from. 

1. Make a docs-related branch:

   ```bash
   git checkout -b docs/pII_open_systems_spinchain
   ```
2. Edit the relevant markdown file in the `docs/` directory. Within the `docs/` directory, we are only concerned with the `docs/qflux/`, which contains separate directories for each module, as well as a global `images/` directory. 

   ```bash
   vim docs/qflux/Open_Systems/spin_chain.md
   ```
3. Preview the documentation to ensure all markdown is correctly rendered within `mkdocs`. You can create a live preview server to see your changes as you edit by running:

   ```bash
   mkdocs serve
   ```
4. If a new page, ensure it is added to the navigation section of `qflux/mkdocs.yml`:
   ```yml
   nav: # Note: All paths must be relative to the docs dir
     - Open Systems:
       - Spin Chain Demo: 'qflux/Open_Systems/spinchainOpen.md'
   ```
5. Commit your changes:

   ```bash
   git add docs/qflux/Open_Systems/spin_chain.md
   git commit -m 'DOCS: Added docs on the spin chain example with Lindblad'
   git push
   ```

> 📌 All documentation should be concise and **not a copy-paste of manuscripts or tutorials**. Use clear markdown formatting. For ease of editing, you can use your favorite markdown editor, IDE, or Jupyter Notebooks to write your initial markdown and then refine/ensure everything renders properly by using your favorite text-editor with the live interactive preview server created via `mkdocs serve`.

---

## Commiting the Update

Follow [conventional commit](https://www.conventionalcommits.org/en/v1.0.0/) guidelines:

* `FEAT:` for new features
* `FIX:` for bug fixes
* `DOCS:` for documentation
* `TEST:` for testing additions
* `REFACTOR:` for code improvements

Example:

```bash
git commit -m "FEAT: Add Lindblad trajectory sampling module"
```

---

## Pull Request

1. Push your feature branch:

   ```bash
   git push origin feat/qumode_dynamics
   ```
2. Open a PR against the `main` branch via GitHub.
3. Provide a meaningful title and description.
4. Link related issues using keywords like `Closes #42`.

---

## Merging a Pull Request

* All PRs should be reviewed and approved by at least one maintainer.
* Ensure the CI (tests + docs build) passes.

---

## Keeping Your Fork Updated

```bash
git remote add upstream https://github.com/batistagroup/qflux.git
git pull upstream main
git push origin main
```

---

Thank you again for contributing!
The QFlux development team appreciates your support.

---
