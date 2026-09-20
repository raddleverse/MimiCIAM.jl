# Baseline Comparison Tests: MimiCIAM dev to MimiCIAM stable
#
# TODO testing against GAMS occurs for now in the GAMS_BaselineComparisons.ipynb
# notebook, the remaining work is to bring those numerical comparisons into this
# repository's scripts.

@testitem "Baseline Comparison: MimiCIAM dev to MimiCIAM stable" setup=[ComparisonData] begin

    using CSV
    using DataFrames

    # provide a directory holding stable julia validation results
    jl_validation_outputdir = joinpath(@__DIR__, "..", "data", "validation_data", "julia")

    files = readdir(jl_validation_outputdir)
    filter!(i -> (i != "desktop.ini" && i != ".DS_Store" && i != "xsc.csv"), files)

    for (i, file) in enumerate(files)

        println("Comparing current CIAM file $(i): $(file) ... to validation version")

        # load data
        expected = CSV.read(joinpath(jl_validation_outputdir, file), DataFrame)
        current = CSV.read(joinpath(comparison_outputdir, file), DataFrame)

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
