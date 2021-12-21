## **Combining Set Propagation with Finite Element Methods for Time Integration in Transient Solid Mechanics Problems**

This repository contains the repeatibility evaluation for the article **Combining Set Propagation with Finite Element Methods for Time Integration in Transient Solid Mechanics Problems** by

- **Marcelo Forets** ([@mforets](http://github.com/mforets)), Departamento de Matemática y Aplicaciones, Centro Universitario Regional del Este, Universidad de la República, Maldonado, Uruguay,
- **Daniel Freire Caporale** ([@dfcaporale](http://github.com/dfcaporale)), Instituto de Fı́sica, Facultad de Ciencias, Universidad de la República, Montevideo, Uruguay), and
- **Jorge M. Pérez Zerpa** ([@jorgepz](http://github.com/jorgepz)), Instituto de Estructuras y Transporte, Facultad de Ingenierı́a, Universidad de la República, Montevideo, Uruguay.

Publication: Forets, Marcelo, Daniel Freire Caporale, and Jorge M. Pérez Zerpa. ["Combining set propagation with finite element methods for time integration in transient solid mechanics problems."](https://www.sciencedirect.com/science/article/abs/pii/S0045794921002212) Computers & Structures 259 (2022): 106699.

The arXiv pre-print can be found here: [arXiv:2105.05841](https://arxiv.org/abs/2105.05841) (2021).

<details>
<summary>Click to read the Abstract.</summary>
  <center>
The Finite Element Method (FEM) is the gold standard for spatial discretization in numerical simulations for a wide spectrum of real-world engineering problems.
	Prototypical areas of interest include linear heat transfer and linear structural dynamics problems modeled with linear partial differential equations (PDEs).
	While different algorithms for direct integration of the equations of motion exist, exploring all feasible behaviors for varying loads, initial states and fluxes in models with large numbers of degrees of freedom remains a challenging task.
	In this article we propose a novel approach, based in set propagation methods and motivated by recent advances in the field of Reachability Analysis.
	Assuming a set of initial states and input states, the proposed method consists in the construction of a union of sets (flowpipe) that enclose the infinite number of solutions of the spatially discretized PDE.
	We present the numerical results obtained in five examples to illustrate the capabilities of the approach, and compare its performance against reference numerical integration methods.
	We conclude that, for problems with single known initial conditions, the proposed method is accurate.
	For problems with uncertain initial conditions included in sets, the proposed method can compute all the solutions of the system more efficiently than numerical integration methods.
    </center>
</details>

## How to run this package

This package requires [Julia](http://julialang.org/). Once installed, run the script `runall.jl` to execute all benchmarks:

```julia
$ julia runall.jl
```
If the execution is successful, a folder `output` will be created, with subfolders
containing the figures for each model.

The Lambs wave propagation problem (problem 3) is not run by default. Pass the flag `LONG` to run that model as well, as in `julia runall.jl LONG`. It may take several hours.

*Note.* For reference, all the figures appearing in the article can be found in the folder [fig](https://github.com/JuliaReach/SetPropagation-FEM-Examples/tree/main/fig).

