using Mimi
using MimiCIAM
using Test
using Query
using RData
using StatsBase
using CSV
using DataFrames
using NetCDF

# load utils functions to write out GAMs and CIAM outputs for baseline comparisons
include("GAMSCIAM_comparison_utils.jl")
include("MimiCIAM_comparison_utils.jl")

@testset begin "Unit Testing"

end

@testset begin "Baseline Comparison: MimiCIAM dev to MimICIAM stable"

end

@testset begin "Baseline Comparison: MimiCIAM dev to GAMS CIAM"

end
