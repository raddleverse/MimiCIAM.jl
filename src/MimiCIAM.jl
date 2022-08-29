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
using Missings

include("slrcost.jl")
include("slrcost_GAMSmatch.jl") # useful for testing against GAMS
include("ciamhelper.jl")
include("ciam.jl")
include("lslr_mapping.jl") # currently setting all lslr = 0 for testing
#include("stub_no_slr.jl")  # a stub all SLR=0 case for testing/hypotheticals

function __init__()
    register(DataDep(
        "BRICK fingerprints",
        """
        Some BRICK fingerprints downloaded from https://github.com/scrim-network/BRICK/raw/master/fingerprints/FINGERPRINTS_SLANGEN_Bakker.nc
        """,
        "https://github.com/scrim-network/BRICK/raw/master/fingerprints/FINGERPRINTS_SLANGEN_Bakker.nc"))
end


end
