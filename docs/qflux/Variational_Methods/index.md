# Variational Methods Module Overview and User Guide

![Logo](../../img/qflux-logo.png)

## Overview

In this section, we outline the main functionality of the `variational_methods` module. 

The `qflux.variational_methods` module aims to provide variational methods for computing dynamics and can be utilized for both open and closed systems that the user is interested in. 

First, we will provide some conceptual explanations that provide the user with a necessary background to understand the code. Then we provide some illustrative examples that demonstrate how the code can be used. Finally, we provide the source code as an API reference to the source code.

## Examples and Introductory Concepts

The first example walks through variational protocols that can be employed for arbitrary real-time evolution `VarQRTE` and imaginary time evolution VarQITE with a simple model system. 

The second example walks through the application of unrestricted adaptive variational quantum dynamics (UAVQD) as applied to a simple instance of the amplitude damping model.

The third example takes things one step further and applies the AVQD protocol in the Stochastic Schrodinger Equation picture as applied to the Fenna-Matthews-Olson (FMO) Complex.

- [Example: Variational Quantum Time Evolution](varQTE.md)
- [Example: Unrestricted Adaptive Variational Quantum Dynamics in an Amplitude Damping Channel](Vectorized_Adaptive.md)
- [Example: Stochastic Schrodinger Equation for Open Systems](trajectory_FMO.md)

