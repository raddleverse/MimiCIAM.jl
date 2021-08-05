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

### Keyword Arguments

The get_model() function has the following signature:
```julia
get_model(;
    initfile = nothing,
    fixed::Bool = false,
    t::Int = 20,
    noRetreat::Bool = false,
    allowMaintain::Bool = true, 
    popinput::Int = 0)
```
which includes several optional keyword arguments to customize the CIAM model you wish to run:
- `initfile` (default = "data/batch/init.csv") : takes a path to a initilization file used to set several parameters (described below) and defaulting to 

- `t` (default = 20): the number of timesteps to run
- `popinput` (default = 0): a socioeconomic parameter that specifies the population data source such with the following options, noting that as of now 1 and 2 are temporarily disabled so 0 is the only option: 0 (default), 1 (Jones & O'Neill, 2016), or 2 (Merkens et al, 2016) 

_we do not recommend altering these without consultation with the authors as changes from the default are experimental_
- `noRetreat` (default = false): a model parameter that specifies if retreat is allowed, such that if the parameter is true, segments will either protect or not adapt, but never retreat.
- `fixed` (default = false): a model parameter that specifies if you want to run the model as fixed (`true`) or flexible (false) with respect to adaptation
- `allowMaintain` (default = true): a model parameter that specifies if maintaining defenses is an option, such that if the parameter is true segments will have the option to maintain current defenses
- 
### Initialization File

The `initfile` parameter above takes a path to a file that must be specially formatted as the `init.txt` file at "data/batch/init.csv":
```
run_name,lslr,subset,ssp,ssp_simplified
base,lsl_rcp85_p50.csv,false,IIASAGDP_SSP5_v9_130219,5
```
This file will indicate the data to import for a given run, the bulk of this work being done in `import_model_data(lslfile, sub, ssp, ssp_simplified, popinput)`.  The file contains several parameters:

- `run_name` (default = base): the name of the run, can be used in labeling and results file production
- `lslr` (default = "lsl_rcp85_p50.csv"): the filename of the file used for lslr settings, which must be available in "data/lslr"
- `subset` (default = false): the list of of segment IDs to run the model for, where false indicates running all segments
- `ssp` (default = "IIASAGDP_SSP5_v9_130219"): the full SSP name that will provide several socioconomic parameters, see the names after "pop" and "ypcc" in "data/ssp" for options
- `ssp_simplified` (default = 5): the integer representing the SSP (1-5)

In order to make creation of such a file easier, we provide an (unexported) file creation function `MimiCIAM.write_init_file(run_name::String, outputdir::String, init_settings::Dict)` which writes the initialization file for a specificied `run_name` into `outputdir` using init_settings
found in `init_settings`.

Note that `init_settings` is a Dictionary with one entry per parameter, best shown through the following example:

```
run_name = "ctrl+SSP5"

init_settings = Dict(
        :init_filename   => string("$run_name", "_init.csv"),
        :lslrfile        => "lsl_rcp85_p50.csv",
        :subset          => false,
        :ssp             => "IIASAGDP_SSP5_v9_130219",
        :ssp_simplified  => 5
    )
    
MimiCIAM.write_init_file(run_name, outputdir, init_settings)
```
