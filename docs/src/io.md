# Input, Output, and Settings

These functions are responsible for reading in all model configurations (passed
by config file or commandline), administrating them during a run, and printing
any output.

## constants.jl

This file defines constants that are needed for calculations involving the
metabolic theory of ecology.

Currently includes `boltz` (the Boltzmann constant) and `act` (the activation
energy).

## defaults.jl

```@autodocs
Modules = [GeMM]
Pages = ["defaults.jl"]
```

## input.jl

```@autodocs
Modules = [GeMM]
Pages = ["input.jl"]
```

## output.jl

```@autodocs
Modules = [GeMM]
Pages = ["output.jl"]
```
