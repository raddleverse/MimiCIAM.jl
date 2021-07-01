using Test
using MimiCIAM

using Query
using RData
using StatsBase
using CSV
using DataFrames
using NetCDF

# write out this working current version of CIAM's output files - note that this can 
# be very slow ... especially the global aggregation ... TODO create a more 
# compact version for quick testing 
include(joinpath(@__DIR__, "write_MimiCIAM_comparison_files.jl"))


@testset begin "Baseline Comparison: MimiCIAM dev to MimICIAM stable"
    
    # a first check here can be against the results found in the folder 
    # test/validation data/MimiCIAM/baselineComparison_07012021
    # which was main branch at that point, then we can do more numerical tests 
    # below

    # TODO: numerical tests against saved MimiCIAM validation data

end

@testset begin "Baseline Comparison: MimiCIAM dev to GAMS CIAM"

    # we use the baselineComparison.ipynb notebook to (1) make graphs and tables 
    # for direct comparison and (2) print out intermediary test files, so first
    # go there for tests, and then more numerical tests will be done below

    # TODO: numerical tests against saved gamsCIAM validation data

end