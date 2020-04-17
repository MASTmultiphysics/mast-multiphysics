# MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit

[![Build Status](https://travis-ci.com/MASTmultiphysics/mast-multiphysics.svg?branch=master)](https://travis-ci.com/MASTmultiphysics/mast-multiphysics)

The Multidisciplinary-design Adaptation and Sensitivity Toolkit (MAST) is a sensitivity-enabled, multiphysics FEA tool developed to support computational design and engineering analysis.

In addition to analysis of complex, nonlinear physics on large built-up models, MAST supports efficient analytical gradient/sensitivity calculation using direct and adjoint methods. As a result, it is well suited for gradient-based multidisciplinary design optimization processes.

MAST is developed at the [Computational Dynamics and Design Laboratory](http://bhatia.ae.msstate.edu) at Mississippi State University in collaboration with the [Air Force Research Laboratory (AFRL)](https://www.afresearchlab.com) Multidisciplinary Science and Technology Center (MSTC).

The MAST website is [https://mastmultiphysics.github.io](https://www.mast-multiphysics.com) and includes examples/tutorials, API documentation, and a growing theory/users guide.

## Capabilities
An abbreviated list of capabilities is given below.
- 1D, 2D, and 3D analysis in multiple disciplines
- Heat Transfer
  - 1D/2D/3D nonlinear heat conduction
  - Radiation and convection boundary conditions
  - Temperature-dependent material properties
- Structures (Elasticity)
  - Beam/plate structural and 2D/3D continuum elements
  - Loading: concentrated force, surface pressure, & thermoelastic
  - Single-point boundary conditions and tie/connector constraints
  - Nonlinear static and transient analysis
  - Modal vibration analysis (including about a nonlinear equilibrium)
  - Bifurcation buckling
  - Advanced continuation and load-stepping in nonlinear analysis
  - Nastran Bulk Data mesh input using [pyNastran](http://pynastran-git.readthedocs.io/en/latest/)
    - other mesh formats also available including Exodus-II and manual definition
- Fluids
  - SU/PG discretization of compressible Euler equations
  - Small-disturbance linearized time-domain and frequency-domain solvers for Euler equations
  - SU/PG discretization of compressible Navier-Stokes equations (experimental)
- Fluid-structure Interaction (FSI)
  - Small-disturbance flutter solution through coupling of structural and fluid discretizations
  - Time-accurate fluid-structure interaction (experimental)
- Aeroelasticity
  - U-g flutter solver with mode tracking
  - Time-domain flutter solver for piston-theory aerodynamics
- Analytical sensitivity analysis for nearly all analysis capabilities
- Interfaces to multiple optimizers (GCMMA, DOT, NPSOL)
- Topology Optimization (TO)
  - Level-set based approaches
  - Density-based approaches

## Installation
Detailed installation/build instructions for MAST and its dependencies are [located on the website](https://www.mast-multiphysics.com/_install.html) or can be located in the source repository in `doc\install\*.dox`.

### Submodules
To keep the size of the main MAST repository smaller, a git submodule is used 
to store large media/assets such as images and animations used for the documentation 
in a separate repo (doc/assets). To build the documentation locally, you must update
the submodule. To do this, simply run the following commands from inside the root
level of this main repository:
```
git submodule init
git submodule update
```

## Copyright
Copyright (C) 2013-2020  Manav Bhatia and MAST authors

## Clearances
Contributions to MAST by AFRL have been cleared for public release by case numbers: 88ABW-2016-5689, 88ABW-2017-5258, 88ABW-2020-0297.