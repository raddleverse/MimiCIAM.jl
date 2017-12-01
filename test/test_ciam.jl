# Catherine Ledna
# November 8, 2017
#
# Test code for CIAM (Coastal Impact and Adaptation Model)
# Adapted from Delavane Diaz (2016) 
#-------------------------------------------------------------------------------
# CIAM computes adaptation costs from sea level rise (slr) in a geographically 
# disaggregated cost-minimization framework. For details and documentation, 
# see Diaz (2016) and https://github.com/delavane/CIAM/.  
#-------------------------------------------------------------------------------

using Mimi
include("../src/ciam.jl")
include("../src/helper.jl")

# For now, evaluating multi segment coastal protection option as proof of concept
# In progress  
# 1. Country to segment, segment to country mapping issue
# 2. Comparison with GAMS version 
#-----------------------------------
# Data Integration
#-----------------------------------

##ISSUE: Mimi does not like arrays of strings
# To get around this, make xsc_ind

# Construction of XSC_IND from XSC: 
# 0. Get XSC from CSV
# TODO functionalize this
datadir = joinpath("..", "data")
parameters = loadparametersciam(datadir)

file="xsc.csv"
xsc_params = Dict{Any, Any}(lowercase(splitext(file)[1]) => readdlm(joinpath(datadir,file), '\r' ))
rgn_seg_params = prepxsc(xsc_params)
xsc_ind_full = rgn_seg_params[1]
allrgns = sort(rgn_seg_params[2])

allsegs = sort(rgn_seg_params[3])
xsc_char_full = rgn_seg_params[4]

include("../src/helper.jl")
testparams = Dict{Any,Any}("countryarea" => parameters["countryarea"], "refpopdens" => parameters["refpopdens"], 
                    "pop" => parameters["pop"], "ypcc" => parameters["ypcc"], "cci" => parameters["cci"],
                    "gtapland" => parameters["gtapland"], "data"=> parameters["data"])

a = parse_ciam_params!(testparams, allrgns, allsegs)
atpers = [1, 4, 8, 12, 16, 20]

# -----------------
# Model Testing
# -----------------
# To do: move non-CSV variables into CSV folder and fully streamline data integration 

m = Model()
setindex(m, :time, 20)
setindex(m, :adaptPers, length(atpers))
setindex(m, :regions, allrgns)
setindex(m, :segments, allsegs)

addcomponent(m, ciam)

# Time params
setparameter(m, :ciam, :ntsteps, 20)       
setparameter(m, :ciam, :tstep, 1.)
setparameter(m, :ciam, :at, atpers) # Testing on t=1 and t=200 to make sure boundaries are right

# Adaptation Parameters
adapt = [1,10,100,1000,10000]
setparameter(m, :ciam, :adaptOptions, adapt)
setparameter(m, :ciam, :fixed, true)

# Socioeconomic params (mostly made up or taken from spreadsheet)
setparameter(m, :ciam, :refpopdens_usa, 50.4)
setparameter(m, :ciam, :ypc_usa, collect(40000.:1000.:60000.))

# Land Params
setparameter(m, :ciam, :landinput, false)
setparameter(m, :ciam, :dvbm, 5.376)
setparameter(m, :ciam, :kgdp, 3.)
setparameter(m, :ciam, :discountrate, 0.04)
setparameter(m, :ciam, :depr, 1.)

# Retreat Params
setparameter(m, :ciam, :mobcapfrac, 0.25)
setparameter(m, :ciam, :movefactor, 1.)
setparameter(m, :ciam, :capmovefactor, .1)
setparameter(m, :ciam, :democost, .05)

# Protection params
setparameter(m, :ciam, :pcfixed, 0.3)
setparameter(m, :ciam, :mc, 0.02)
setparameter(m, :ciam, :pc0, 6.02)

# Storm Params
setparameter(m, :ciam, :floodmortality, 0.01)

# Wetland params
setparameter(m, :ciam, :wbvm, 0.376)
setparameter(m, :ciam, :wmaxrate, 0.01)

setparameter(m, :ciam, :xsc, xsc_ind_full)

# Dummy lslr variable
lslr_allsegs = transpose(collect(1.5:.1:3.4))
lslr_allsegs = repeat(lslr_allsegs, outer = (12148,1))

setparameter(m, :ciam, :lslr, lslr_allsegs)

tp = Dict{Any,Any}("pop" => testparams["pop"],
"refpopdens"=> transpose(testparams["refpopdens"]))
setleftoverparameters(m, testparams)

include("../src/ciam.jl")

@time run(m)



