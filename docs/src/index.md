```@meta
CurrentModule = GeMM
```

# Introduction

This is the documentation for the Island Speciation Model, a.k.a. Genetically
explicit Metacommunity Model (GeMM). *(Yeah, we're still working on the name...)*

The aim of this model is to create a virtual island ecosystem that can be used to
explore ecological and evolutionary hypotheses *in silico*. It is genetically
and spatially explicit, with discrete space and time.

This documentation is generated from the source code using Julia's inbuilt
`Documenter` module. It is sorted first by functionality, second by source code
file.

## Running the model

The functions in the `run_simulation.jl` file are used to start a simulation run:

```@autodocs
Modules = [GeMM]
Pages = ["run_simulation.jl"]
```

`rungemmparallel.jl` bundles these for quick access. From the commandline, call:

```
> julia -p <cores> rungemmparallel.jl -c <configs>
```

where `<cores>` is the number of processors you want to make available to Julia
(one processor per simulation max), and `<configs>` is a comma-separated list of
configuration files that are to be processed.

*Last updated: 2019-02-15 (commit 27a93d0)*  
