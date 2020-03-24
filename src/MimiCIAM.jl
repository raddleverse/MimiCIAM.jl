module MimiCIAM

using Mimi
using DelimitedFiles
using DataFrames
using Query
using CSV
using Dates
using NetCDF
using Statistics
using StatsBase

include("slrcost.jl")
include("ciamhelper.jl")
include("ciam.jl")
include("brickLSL.jl")

end
