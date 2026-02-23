# GridapRobustDarcyForchheimer

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ruizbaier.github.io/GridapRobustDarcyForchheimer.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ruizbaier.github.io/GridapRobustDarcyForchheimer.jl/dev/)
[![Build Status](https://github.com/ruizbaier/GridapRobustDarcyForchheimer.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ruizbaier/GridapRobustDarcyForchheimer.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ruizbaier/GridapRobustDarcyForchheimer.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ruizbaier/GridapRobustDarcyForchheimer.jl)

## How to cite

This repository accompanies the paper [[1]],  submitted to SIAM Journal on Numerical Analysis. 

```
@article{dhpr_DarcyForchheimer26,
    author = {Das, Rishi and Hutridurga, Harsha and Pani, Amiya K. and Ruiz-Baier, Ricardo},
    doi = {10.48550/arXiv.2510.24527},
    journal = {arXiv preprint 2510.24527},
    year    = {2026},
    title   = {Robust stability and preconditioning of Darcy--Forchheimer equations},
}
```

It contains the implementation of the mixed method and preconditioners proposed in the paper, as well as the tests that verify convergence and robustness with respect to the physical and discretisation parameters.

The implementation relies on the [Gridap](https://gridap.github.io/Gridap.jl/stable/) finite element library. Visualisation is done with [Paraview](https://paraview.org). For the mesh generation and manipulation we use the library [GMSH](https://gmsh.info). 

