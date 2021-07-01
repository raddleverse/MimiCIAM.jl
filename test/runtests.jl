using Mimi
using MimiCIAM
using Test
using Query
using RData
using StatsBase
using CSV
using DataFrames
using NetCDF

write_out_GAMS = false # run only once to produce validation files
write_out_CIAM = true

# write out GAMs and CIAM outputs for baseline comparisons
write_out_GAMS && include(joinpath(@__DIR__, "write_GAMSCIAM_comparison_files.jl"))
write_out_CIAM && include(joinpath(@__DIR__, "write_MimiCIAM_comparison_files.jl"))

@testset begin "Unit Testing"

end

@testset begin "Baseline Comparison: MimiCIAM dev to MimICIAM stable"

end

@testset begin "Baseline Comparison: MimiCIAM dev to GAMS CIAM"

end
