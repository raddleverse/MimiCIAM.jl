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
include("../src/ciam_rgn.jl")

# For now, evaluating single segment coastal protection option as proof of concept
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
setparameter(m, :ciam, :length, 86.08)
setparameter(m, :ciam, :depr, 1.)

# Retreat Params
setparameter(m, :ciam, :mobcapfrac, 0.25)
setparameter(m, :ciam, :movefactor, 1.)
setparameter(m, :ciam, :capmovefactor, .1)
setparameter(m, :ciam, :democost, .05)


# Protection params
setparameter(m, :ciam, :pc, 2.)
setparameter(m, :ciam, :pcfixed, 0.3)
setparameter(m, :ciam, :mc, 0.02)
setparameter(m, :ciam, :pc0, 6.02)

# Surge exposure params
setparameter(m, :ciam, :pσ₀, 2.05)
setparameter(m, :ciam, :pσ₀coef, 26.4)
setparameter(m, :ciam, :pσ₁, 0.06)
setparameter(m, :ciam, :pσ₂, 18.5)
setparameter(m, :ciam, :rσ₀, [26.4 26.4 26.4 26.4])
setparameter(m, :ciam, :rσ₁, [0.06 .06 .06 .06])
setparameter(m, :ciam, :rσ₂, [18.5 18.5 18.5 18.5])


setparameter(m, :ciam, :floodmortality, 0.01)

# Wetland params
setparameter(m, :ciam, :wbvm, 0.376)
setparameter(m, :ciam, :wetlandarea, [30. 30. 30. 30.])
setparameter(m, :ciam, :wmaxrate, 0.01)

# Coast area params
areaparams = [152., 4., 2., 1., 1., 1.333, 1.333, 1.333, 1., 1., 1., 1., 1.5, 1.5, 1.5]
areaparams_all = [transpose(areaparams); transpose(areaparams); transpose(areaparams); transpose(areaparams)]
setparameter(m, :ciam, :areaparams, areaparams_all)

# Slr params
setparameter(m, :ciam, :lslr, lslr)
setparameter(m, :ciam, :surgeExposure, [0., 0.1, 0.2, 0.3, 3.3] )
setparameter(m, :ciam, :xsc, xsc_ind)

@time run(m)


# Sanity Checks
m[:ciam, :AdaptationDecision]
m[:ciam, :AdaptationLevel]
m[:ciam, :AdaptationCost]

m[:ciam, :ProtectCost]
m[:ciam, :RetreatCost]
m[:ciam, :NoAdaptCost]

### TO DO THIS WEEK
# 1. Protect - last period, what's up?---done
# 2. Sanity checks - protect, retreat, no adapt---done
#       2a. Decision - protectlevel, retreat level params?---done
#       2b. Fixed vs flexible
# 3. Thorough debug - track calculation--- today
#   3a. Upload to google drive today - done
# 4. Timing - single segment
# 5. Multisegment integration 
#   5a. Github upload 
# 6. (Simultaneously) data integration
# 7. Comparison with Delavane


### Case 2: Same params, but different AT boundary cases
m2 = Model()
setindex(m2, :time, 20)

addcomponent(m2, ciam)

# Time params
setparameter(m2, :ciam, :tstep, 1.)
setparameter(m2, :ciam, :at, [3, 5, 10, 18]) # Testing on t=3 and t=18 to make sure boundaries are right
setparameter(m2, :ciam, :ntsteps, 20.)       

# Socioeconomic params (mostly made up or taken from spreadsheet)
setparameter(m2, :ciam, :pop_country, collect(350.:10.:550.))
setparameter(m2, :ciam, :refpopdens_country, 1.8)
setparameter(m2, :ciam, :refpopdens_usa, 50.4)

setparameter(m2, :ciam, :popdens1_seg, 1.8)
setparameter(m2, :ciam, :ypc_country, collect(35000.:1000.:55000.))
setparameter(m2, :ciam, :ypc_usa, collect(40000.:1000.:60000.))

# Land Params
setparameter(m2, :ciam, :landinput, 0.)
setparameter(m2, :ciam, :gtapland, 0.13)
setparameter(m2, :ciam, :dvbm, 5.376)
setparameter(m2, :ciam, :kgdp, 3.)
setparameter(m2, :ciam, :discountrate, 0.04)
setparameter(m2, :ciam, :length, 86.08)

# Protection params
setparameter(m2, :ciam, :pc, 2.)
setparameter(m2, :ciam, :pcfixed, 0.3)
setparameter(m2, :ciam, :mc, 0.02)
setparameter(m2, :ciam, :pc0, 6.02)

# Surge exposure params
setparameter(m2, :ciam, :pσ₀, 2.05)
setparameter(m2, :ciam, :pσ₀coef, 26.4)
setparameter(m2, :ciam, :pσ₁, 0.06)
setparameter(m2, :ciam, :pσ₂, 18.5)

setparameter(m2, :ciam, :floodmortality, 0.01)

# Wetland params
setparameter(m2, :ciam, :wbvm, 0.376)
setparameter(m2, :ciam, :wetlandarea, 30.)

# Slr params
setparameter(m2, :ciam, :lslr, lslr)
setparameter(m2, :ciam, :slr10, 0.1)
setparameter(m2, :ciam, :slr100, 0.2)
setparameter(m2, :ciam, :slr1000, 0.3)
setparameter(m2, :ciam, :slr10000, 3.3)


# Sanity check for model outcomes
m2[:ciam, :H10]
m2[:ciam, :H100]
m2[:ciam, :H1000]
m2[:ciam, :H10000]

m2[:ciam, :ProtectCost]
m2[:ciam, :ProtectLevel]
m2[:ciam, :AdaptationOption]
m2[:ciam, :AdaptationLevel]
m2[:ciam, :AdaptationCost]


# Evaluate Model Results/Outcomes
runtimes = ones(1.:100.)
runtimes2 = ones(1.:100.)
for i in collect(1:100)
    runtimes[i] = (@timed run(m))[2]
    runtimes2[i] = (@timed run(m2))[2]
end

mean(runtimes)
mean(runtimes2)