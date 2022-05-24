using Test
using CSV
using DataFrames
using MimiCIAM

include("test_utils.jl")

##==============================================================================
## Unit Testing

@testset "Unit Testing" begin
    m = MimiCIAM.get_model()
    run(m)
end

##==============================================================================
## Baseline Comparison Tests - Gather Data (takes ~20 minutes)

# provide a directory to save julia results
jl_outputdir = joinpath(@__DIR__, "..", "output", "BaselineComparisonTests")

# write out the current results using the same random subset of segments used
# for the validation data
write_MimiCIAM_comparison_files(jl_outputdir, subset = "random10.csv")

##==============================================================================
## Baseline Comparison Tests: MimiCIAM dev to MimICIAM stable

@testset "Baseline Comparison: MimiCIAM dev to MimICIAM stable" begin

    # provide a directory holding stable julia validation results
    jl_validation_outputdir = joinpath(@__DIR__, "..", "data", "validation_data", "julia")

    files = readdir(jl_validation_outputdir)
    filter!(i->(i!="desktop.ini" && i!=".DS_Store" && i!="xsc.csv"), files)

    for (i, file) in enumerate(files)

        println("Comparing current CIAM file $(i): $(file) ... to validation version")

        # load data
        expected = CSV.read(joinpath(jl_validation_outputdir, file), DataFrame)
        current = CSV.read(joinpath(jl_outputdir, file), DataFrame)

        # sort data
        fields = DataFrames.names(expected)
        value_field = Symbol.(fields[end])
        sort_fields = Symbol.(fields[1:end-1])

        sort!(expected, sort_fields)
        sort!(current, sort_fields)

        # compare values
        diffs = abs.(expected[!, value_field] .- current[!, value_field])
        @test maximum(diffs) <= 1e-8 # handles rounding errors from different Excel versions
    end
end

##==============================================================================
## Baseline Comparison Tests: MimiCIAM dev to GAMS

@testset "Baseline Comparison: MimiCIAM dev to GAMS CIAM" begin

    # testing against GAMS ocurrs for now in the GAMS_BaselineComparisons.ipynb
    # notebook, TODO is to bring the numerical comparisons into this repository's
    # scripts
end
