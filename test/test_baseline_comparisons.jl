using Test
using MimiCIAM

using Query
using RData
using StatsBase
using CSV
using DataFrames
using NetCDF

# TODO - avoid all the local path nonsense :)

##==============================================================================
## Gather Data (~30 minutes)

# provide a directory to save julia results, defaults to output folder but we
# recommend providing a local path
# jl_outputdir = joinpath(@__DIR__, "..", "output", "results-jl")
jl_outputdir = "/Users/lisarennels/JuliaProjects/CIAMPaper/local-data/jl-outputs-curr/raw"

write_MimiCIAM_comparison_files(outputdir = jl_outputdir)

##==============================================================================
## Test (~30 minutes)

@testset "Baseline Comparison: MimiCIAM dev to GAMS CIAM" begin

    # STEP 1
    # The first main test is to use the baselineComparison.ipynb notebook, which 
    # will take in directory names and do comparisons

    # STEP 2
    # TODO: numerical tests against saved gamsCIAM validation data

    # provide a directory holding julia validation results, defaults to validation 
    # folder (currently empty)
    # gams_validation_outputdir = joinpath(@__DIR__, "validation_data", "gams")
    gams_validation_outputdir = "/Users/lisarennels/JuliaProjects/CIAMPaper/local-data/gams-outputs"

end

# For now we do not have exact numerical tests
@testset "Baseline Comparison: MimiCIAM dev to MimICIAM stable" begin

    # provide a directory holding gams validation results, defaults to validation 
    # folder (currently empty)
    # jl_validation_outputdir = joinpath(@__DIR__, "validation_data", "julia")
    jl_validation_outputdir = "/Users/lisarennels/JuliaProjects/CIAMPaper/local-data/jl-outputs-05042021"

    files = [
                "ctrl+noConstrFix_global_85p50ssp0fixed.csv", 
                "ctrl+noConstrFix_seg_85p50ssp0fixed_optimal.csv", 
                "ctrl+noConstrFix_seg_85p50ssp0fixed.csv"
            ]

    for (i, file) in enumerate(files)

        label_fields = (file == 1) ? [:time, :regions, :segments, :variable, :level] : [:time, :variable, :level]
        value_field = (file < 3) ? :value : :OptimalCost

        # load validation data
        expected = CSV.read(joinpath(jl_validation_outputdir, file), DataFrame)
        sort!(expected, label_fields)

        # load current data
        current = CSV.read(joinpath(jl_outputdir, file), DataFrame)
        sort!(expected, label_fields)

        # check the organization
        for col in label_fields
            @test filter(x -> !ismissing(x), expected[!, col]) == filter(x -> !ismissing(x), current[!, col])
        end

        # compare values
        diffs = abs.(expected[!, value_field] .- current[!, :value_field])
        @test maximum(diffs) <= 1e-9

    end

end

