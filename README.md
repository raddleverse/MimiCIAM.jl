# MimiCIAM.jl

This is a work-in-progress respository for a Julia-Mimi implementation the CIAM model adapted from Diaz, 2016.

## Preparing the Software Environment

Your first step is to install MimiCIAM.jl itself, and to do so you need to run the following command at the julia package REPL:

```julia
pkg> add https://github.com/BRICK-SLR/MimiCIAM.jl.git
```

You probably also want to install the Mimi package into your julia environment, so that you can use some of the tools in there:

```julia
pkg> add Mimi
```
## Running the Model

The model uses the Mimi framework and it is highly recommended to read the Mimi  documentation first to understand the code structure. The basic way to access a copy of the default MimiFAIRv2 model and explore the resuts is the following:

The basic way to access a copy of the default MimiCIAM model is the following:

```julia
using MimiCIAM
m = MimiCIAM.get_model()
run(m)
```

The get_model() function currently has the following keyword arguments:
- `initfile`: (default = nothing)
- `fixed`: (default = false)
- `noRetreat`: (default = false)
- `allowMaintain`: (default = true)
- `popinput`: (default = 0)
