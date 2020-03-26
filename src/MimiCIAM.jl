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
using DataDeps

include("slrcost.jl")
include("ciamhelper.jl")
include("ciam.jl")
include("brickLSL.jl")

function __init__()
    register(DataDep(
        "BRICK fingerprints",
        """
        Some BRICK fingerprints downloaded from https://github.com/scrim-network/BRICK/raw/master/fingerprints/FINGERPRINTS_SLANGEN_Bakker.nc
        """,
        "https://github.com/scrim-network/BRICK/raw/master/fingerprints/FINGERPRINTS_SLANGEN_Bakker.nc"))    
end


end
