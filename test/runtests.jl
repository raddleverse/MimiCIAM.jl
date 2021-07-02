using MimiCIAM
using Test
using CSV

# TODO - avoid all the local path nonsense and get things uploaded somewhere :)

include("utils.jl")

##==============================================================================
## Unit Testing

@testset "Unit Testing" begin 
    m = MimiCIAM.get_model()
    run(m)
end

##==============================================================================
## Baseline Comparison Tests - Gather Data (takes ~20 minutes)

# provide a directory to save julia results
# jl_outputdir = joinpath(@__DIR__, "..", "output", "results-jl")
jl_outputdir = "/Users/lisarennels/JuliaProjects/CIAMPaper/local-data/jl-outputs-curr/raw"
write_MimiCIAM_comparison_files(jl_outputdir)

##==============================================================================
## Baseline Comparison Tests: MimiCIAM dev to MimICIAM stable

@testset "Baseline Comparison: MimiCIAM dev to MimICIAM stable" begin

    # provide a directory holding stable julia validation results
    # jl_validation_outputdir = joinpath(@__DIR__, "validation_data", "julia")
    jl_validation_outputdir = "/Users/lisarennels/JuliaProjects/CIAMPaper/local-data/jl-outputs-07012021"

    files = [
                "ctrl+noConstrFix_global_85p50ssp0fixed.csv", 
                "ctrl+noConstrFix_seg_85p50ssp0fixed_optimal.csv", 
                "ctrl+noConstrFix_seg_85p50ssp0fixed.csv"
            ]

    for (i, file) in enumerate(files)

        println("Comparing current CIAM file $(i): $(file) ... to validation version")
        label_fields = (i > 1) ? [:time, :regions, :segments, :variable, :level] : [:time, :variable, :level]
        value_field = (i == 2) ? :OptimalCost : :value

        # load validation data
        expected = CSV.read(joinpath(jl_validation_outputdir, file), DataFrame)

        # load current data
        current = CSV.read(joinpath(jl_outputdir, file), DataFrame)

        # align the dataframe ordering
        sort!(expected, label_fields)
        sort!(current, label_fields)
        for col in label_fields
            @test filter(x -> !ismissing(x), expected[!, col]) == filter(x -> !ismissing(x), current[!, col])
        end

        # compare values
        diffs = abs.(expected[!, value_field] .- current[!, value_field])
        @test maximum(diffs) <= 1e-9
    end
end

##==============================================================================
## Baseline Comparison Tests: MimiCIAM dev to GAMS

@testset "Baseline Comparison: MimiCIAM dev to GAMS CIAM" begin

    # STEP 1
    # The first main test is to use the baselineComparison.ipynb notebook, which 
    # will take in directory names and do comparisons

    # STEP 2
    # TODO: numerical tests against saved gamsCIAM validation data

    # provide a directory holding julia validation results
    # gams_validation_outputdir = joinpath(@__DIR__, "validation_data", "gams")
    gams_validation_outputdir = "/Users/lisarennels/JuliaProjects/CIAMPaper/local-data/gams-outputs"

end
