# Catherine Ledna
# October 6, 2017
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

# For now, evaluating multi segment coastal protection option as proof of concept
# Still not implemented: 
# 1. Country to segment, segment to country mappings
# 2. Data integration/formatting
# 3. Comparison with GAMS version (this is just a proof of concept, coding review portion)


# Dummy lslr variable
lslr = [transpose(collect(1.5:.1:3.4)); transpose(collect(1.4:.1:3.3));
        transpose(collect(0.1:.1:2.0));transpose(collect(0.:.1:1.9))]
adapt = [1,10,100,1000,10000]
atpers = [1, 4, 8, 12, 16, 20]
rgns = ["USA", "CAN"]
segs = ["UnitedStates4670", "UnitedStates4648", "Canada5413", "Canada5287"]

##ISSUE: Mimi does not like arrays of strings
# To get around this, make xsc_ind
#   Will need to automate in future

xsc_ind = Dict(1 => 1, 2 => 1, 3 =>2, 4 =>2)


m = Model()
setindex(m, :time, 20)
setindex(m, :level, length(adapt))
setindex(m, :adaptPers, length(atpers))
setindex(m, :regions, rgns)
setindex(m, :segments, segs)

addcomponent(m, ciam)

# Region / segment params: 
#setparameter(m, :ciam, :rgn_names, rgns)
#setparameter(m, :ciam, :seg_names, segs)


# Time params
setparameter(m, :ciam, :tstep, 1.)
setparameter(m, :ciam, :at, atpers) # Testing on t=1 and t=200 to make sure boundaries are right
setparameter(m, :ciam, :ntsteps, 20)       

# Adaptation Parameters
setparameter(m, :ciam, :adaptOptions, adapt)

setparameter(m, :ciam, :fixed, true)

# Socioeconomic params (mostly made up or taken from spreadsheet)
setparameter(m, :ciam, :pop_country, [transpose(collect(350.:10.:550.)); transpose(collect(250.:10.:450.))])
setparameter(m, :ciam, :refpopdens_country, [1.8; 1.2])
setparameter(m, :ciam, :refpopdens_usa, 50.4)

setparameter(m, :ciam, :popdens1_seg, [1.8;1.5;1.7;.9])
setparameter(m, :ciam, :ypc_country, [transpose(collect(35000.:1000.:55000.)); transpose(collect(25000.:1000.:45000.))])
setparameter(m, :ciam, :ypc_usa, collect(40000.:1000.:60000.))

# Land Params
setparameter(m, :ciam, :landinput, false)
setparameter(m, :ciam, :gtapland, [0.13; 0.12])
setparameter(m, :ciam, :dvbm, 5.376)
setparameter(m, :ciam, :kgdp, 3.)
setparameter(m, :ciam, :discountrate, 0.04)
setparameter(m, :ciam, :length, [86.08; 80.; 80.; 80.])
setparameter(m, :ciam, :depr, 1.)

# Retreat Params
setparameter(m, :ciam, :mobcapfrac, 0.25)
setparameter(m, :ciam, :movefactor, 1.)
setparameter(m, :ciam, :capmovefactor, .1)
setparameter(m, :ciam, :democost, .05)


# Protection params
setparameter(m, :ciam, :pc, [2.; 2.; 2.; 2.])
setparameter(m, :ciam, :pcfixed, 0.3)
setparameter(m, :ciam, :mc, 0.02)
setparameter(m, :ciam, :pc0, 6.02)

# Surge exposure params
setparameter(m, :ciam, :pσ₀, [2.05; 2.05; 2.05; 2.05])
setparameter(m, :ciam, :pσ₀coef, [26.4;26.4;26.4;26.4])
setparameter(m, :ciam, :pσ₁, [0.06; 0.06; 0.06; 0.06])
setparameter(m, :ciam, :pσ₂, [18.5; 18.5; 18.5; 18.5])
setparameter(m, :ciam, :rσ₀, [26.4, 26.4, 26.4, 26.4])
setparameter(m, :ciam, :rσ₁, [0.06, .06, .06, .06])
setparameter(m, :ciam, :rσ₂, [18.5, 18.5, 18.5, 18.5])


setparameter(m, :ciam, :floodmortality, 0.01)

# Wetland params
setparameter(m, :ciam, :wbvm, 0.376)
setparameter(m, :ciam, :wetlandarea, [30.; 30.; 30.; 30.])
setparameter(m, :ciam, :wmaxrate, 0.01)

# Coast area params
areaparams = [152., 4., 2., 1., 1., 1.333, 1.333, 1.333, 1., 1., 1., 1., 1.5, 1.5, 1.5]
areaparams_all = [transpose(areaparams); transpose(areaparams); transpose(areaparams); transpose(areaparams)]
setparameter(m, :ciam, :areaparams, areaparams_all)

# Slr params
setparameter(m, :ciam, :lslr, lslr)
surgeExposure = [0., 0.1, 0.2, 0.3, 3.3]
surgeExposureAll = [transpose(surgeExposure); transpose(surgeExposure); transpose(surgeExposure); transpose(surgeExposure)]

setparameter(m, :ciam, :surgeExposure, surgeExposureAll )
setparameter(m, :ciam, :xsc, xsc_ind)

@time run(m)


# Sanity Checks
m[:ciam, :AdaptationDecision]
m[:ciam, :AdaptationLevel]
m[:ciam, :AdaptationCost]

m[:ciam, :ProtectCost]
m[:ciam, :RetreatCost]
m[:ciam, :NoAdaptCost]
m[:ciam, :R]
m[:ciam, :H]

# Evaluate Model Results/Outcomes
runtimes = ones(1.:100.)
for i in collect(1:100)
    runtimes[i] = (@timed run(m))[2]

end

mean(runtimes)
