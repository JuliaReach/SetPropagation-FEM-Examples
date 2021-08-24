## **Combining Set Propagation with Finite Element Methods for Time Integration in Transient Solid Mechanics Problems** (Repeatability Evaluation package)

This repository contains the repeatibility evaluation for the article **Combining Set Propagation with Finite Element Methods for Time Integration in Transient Solid Mechanics Problems** by **Marcelo Forets** ([@mforets](http://github.com/mforets) from Departamento de Matemática y Aplicaciones, Centro Universitario Regional del Este, Universidad de la República, Maldonado, Uruguay, **Daniel Freire Caporale** ([@dfcaporale](http://github.com/dfcaporale)) from Instituto de Fı́sica, Facultad de Ciencias, Universidad de la República, Montevideo, Uruguay), and **Jorge M. Pérez Zerpa** ([@jorgepz](http://github.com/jorgepz)) from Instituto de Estructuras y Transporte, Facultad de Ingenierı́a, Universidad de la República, Montevideo, Uruguay.

The arXiv pre-print can be found here: [arXiv:2105.05841](https://arxiv.org/abs/2105.05841) (2021).

**Abstract.** *The Finite Element Method (FEM) is the gold standard for spatial discretization in
numerical simulations for a wide spectrum of real-world engineering problems.
Prototypical areas of interest include linear heat transfer and linear structural
dynamics problems modeled with linear partial differential equations (PDEs).
While different algorithms for direct integration of the equations of motion exist,
exploring all feasible behaviors for varying loads, initial states and fluxes in models
with large numbers of degrees of freedom remains a challenging task. In this article
we propose a novel approach, based in set propagation methods and motivated by recent
advances in the field of Reachability Analysis. Assuming a set of initial states and
input states, the proposed method consists in the construction of a union of sets
(flowpipe) that enclose the infinite number of solutions of the spatially discretized PDE.
We present the numerical results obtained in four examples to illustrate the capabilities
of the approach, and draw some comparisons with respect to reference numerical integration
methods. We conclude that the proposed method presents specific and promising advantages,
but the full potential of reachability analysis in solid mechanics problems is yet to be explored.*

## How to run this package

Run the script `runall.jl` to execute all benchmarks:

```julia
$ julia runall.jl
```
If the execution is successful, a folder `output` will be created, with subfolders
containing the figures for each model.

All the figures appearing in the article can be found in the folder [fig](https://github.com/JuliaReach/SetPropagation-FEM-Examples/tree/main/fig).
