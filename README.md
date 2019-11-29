# Accelerating MOEA/D by Nelder-Mead method

This code is for the paper [Accelerating MOEA/D by Nelder-Mead method](https://ieeexplore.ieee.org/document/7969414/)

## Abstract:

The multiobjective evolutionary algorithm based on decomposition (MOEA/D) converts a multiobjective optimization problem into a set of single-objective subproblems, and tackles them simultaneously. In MOEA/D, the offspring generation is a crucial part to increase the convergence of the algorithm and maintain the diversity of the solution set. Currently, the majority of reproduction operators consider the quality of neighborhood exploration, i.e., the capability to distribute along the population structure, while few operators have good capability for subproblem exploitation, i.e., the ability to push solutions forward along the subproblems. To address this issue in this paper, we introduce one of the derivative-free optimization methods, Nelder-Mead simplex (NMS) method, to MOEA/D to accelerate the algorithm convergence. The NMS operator is combined with a differential evolution (DE) operator in the offspring generation. The comparison study demonstrates that calling the NMS operator occasionally can help to accelerate the convergence.

## Code

The MOEA/D framework is implemented in MATLAB, while the Nelder-Mead simplex is implemented in C++. 
So you need "mex" to compile the code in MATLAB as below:

> mex --setup C++

> mex NSS.cpp
