using Test
using MimiCIAM

using Query
using RData
using StatsBase
using CSV
using DataFrames
using NetCDF

##==============================================================================
## Setup

# TODO - add subsets of the data to the validation_data folder so we can avoid
# storing data locally

# provide a directory to save julia results, defaults to output folder but we
# recommend providing a local path

# recommend providing a local `outputdir` to use, but it will default to the output 
# directory below
jl_outputdir = joinpath(@__DIR__, "..", "output", "results-jl")

# provide a directory holding julia validation results, defaults to validation 
# folder (currently empty)
jl_validation_outputdir = joinpath(@__DIR__, "validation_data", "julia")

# provide a directory holding gams validation results, defaults to validation 
# folder (currently empty)
gams_validation_outputdir = joinpath(@__DIR__, "validation_data", "gams")

##==============================================================================
## Gather Data (~30 minutes)

write_MimiCIAM_comparison_files(outputdir = jl_outputdir)

##==============================================================================
## Test (~30 minutes)

@testset "Baseline Comparison: MimiCIAM dev to GAMS CIAM" begin

    # The first main test is to use the baselineComparison.ipynb notebook, which 
    # will take in directory names and do comparisons

    # TODO: numerical tests against saved gamsCIAM validation data

end

# For now we do not have exact numerical tests
@testset "Baseline Comparison: MimiCIAM dev to MimICIAM stable" begin

    # TODO: numerical tests against saved MimiCIAM validation data

end

