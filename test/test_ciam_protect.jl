# Catherine Ledna
# October 6, 2017
#
# Test code for CIAM (Coastal Impact and Adaptation Model)
# Adapted from Delavane Diaz (2016) 
#-------------------------------------------------------------------------------
# CIAM computes adaptation costs from sea level rise (slr) in a geographically 
# disaggregated cost-minimization framework. For details and documentation, 
# see Diaz (2016) and https://github.com/delavane/CIAM/.  
#
# This is a truncated version for testing purposes that only evaluates the coastal
# protection option. 
#-------------------------------------------------------------------------------

using Mimi
include("../src/ciam_protect.jl")

# For now, evaluating single segment coastal protection option as proof of concept
# Still not implemented: 
# 1. Country to segment, segment to country mappings
# 2. Data integration/formatting
# 3. Comparison with GAMS version (this is just a proof of concept, coding review portion)


# Dummy lslr variable
lslr = collect(1.5:.1:3.4)

m = Model()
setindex(m, :time, 20)

addcomponent(m, ciam_protect)

# Time params
setparameter(m, :ciam_protect, :tstep, 1.)
setparameter(m, :ciam_protect, :at, [1, 5, 10, 20]) # Testing on t=1 and t=20 to make sure boundaries are right
setparameter(m, :ciam_protect, :ntsteps, 20.)       

# Socioeconomic params (mostly made up or taken from spreadsheet)
setparameter(m, :ciam_protect, :pop_country, collect(350.:10.:550.))
setparameter(m, :ciam_protect, :refpopdens_country, 1.8)
setparameter(m, :ciam_protect, :refpopdens_usa, 50.4)

setparameter(m, :ciam_protect, :popdens1_seg, 1.8)
setparameter(m, :ciam_protect, :ypc_country, collect(35000.:1000.:55000.))
setparameter(m, :ciam_protect, :ypc_usa, collect(40000.:1000.:60000.))

# Land Params
setparameter(m, :ciam_protect, :landinput, 0.)
setparameter(m, :ciam_protect, :gtapland, 0.13)
setparameter(m, :ciam_protect, :dvbm, 5.376)
setparameter(m, :ciam_protect, :kgdp, 3.)
setparameter(m, :ciam_protect, :discountrate, 0.04)
setparameter(m, :ciam_protect, :length, 86.08)

# Protection params
setparameter(m, :ciam_protect, :pc, 2.)
setparameter(m, :ciam_protect, :pcfixed, 0.3)
setparameter(m, :ciam_protect, :mc, 0.02)
setparameter(m, :ciam_protect, :pc0, 6.02)

# Surge exposure params
setparameter(m, :ciam_protect, :pσ₀, 2.05)
setparameter(m, :ciam_protect, :pσ₀coef, 26.4)
setparameter(m, :ciam_protect, :pσ₁, 0.06)
setparameter(m, :ciam_protect, :pσ₂, 18.5)

setparameter(m, :ciam_protect, :floodmortality, 0.01)

# Wetland params
setparameter(m, :ciam_protect, :wbvm, 0.376)
setparameter(m, :ciam_protect, :wetlandarea, 30.)

# Slr params
setparameter(m, :ciam_protect, :lslr, lslr)
setparameter(m, :ciam_protect, :slr10, 0.1)
setparameter(m, :ciam_protect, :slr100, 0.2)
setparameter(m, :ciam_protect, :slr1000, 0.3)
setparameter(m, :ciam_protect, :slr10000, 3.3)


@time run(m)@


# Sanity check for model outcomes
m[:ciam_protect, :H10]
m[:ciam_protect, :H100]
m[:ciam_protect, :H1000]
m[:ciam_protect, :H10000]

m[:ciam_protect, :ProtectCost]
m[:ciam_protect, :ProtectLevel]
m[:ciam_protect, :AdaptationOption]
m[:ciam_protect, :AdaptationLevel]
m[:ciam_protect, :AdaptationCost]


### Case 2: Same params, but different AT boundary cases
m2 = Model()
setindex(m2, :time, 20)

addcomponent(m2, ciam_protect)

# Time params
setparameter(m2, :ciam_protect, :tstep, 1.)
setparameter(m2, :ciam_protect, :at, [3, 5, 10, 18]) # Testing on t=3 and t=18 to make sure boundaries are right
setparameter(m2, :ciam_protect, :ntsteps, 20.)       

# Socioeconomic params (mostly made up or taken from spreadsheet)
setparameter(m2, :ciam_protect, :pop_country, collect(350.:10.:550.))
setparameter(m2, :ciam_protect, :refpopdens_country, 1.8)
setparameter(m2, :ciam_protect, :refpopdens_usa, 50.4)

setparameter(m2, :ciam_protect, :popdens1_seg, 1.8)
setparameter(m2, :ciam_protect, :ypc_country, collect(35000.:1000.:55000.))
setparameter(m2, :ciam_protect, :ypc_usa, collect(40000.:1000.:60000.))

# Land Params
setparameter(m2, :ciam_protect, :landinput, 0.)
setparameter(m2, :ciam_protect, :gtapland, 0.13)
setparameter(m2, :ciam_protect, :dvbm, 5.376)
setparameter(m2, :ciam_protect, :kgdp, 3.)
setparameter(m2, :ciam_protect, :discountrate, 0.04)
setparameter(m2, :ciam_protect, :length, 86.08)

# Protection params
setparameter(m2, :ciam_protect, :pc, 2.)
setparameter(m2, :ciam_protect, :pcfixed, 0.3)
setparameter(m2, :ciam_protect, :mc, 0.02)
setparameter(m2, :ciam_protect, :pc0, 6.02)

# Surge exposure params
setparameter(m2, :ciam_protect, :pσ₀, 2.05)
setparameter(m2, :ciam_protect, :pσ₀coef, 26.4)
setparameter(m2, :ciam_protect, :pσ₁, 0.06)
setparameter(m2, :ciam_protect, :pσ₂, 18.5)

setparameter(m2, :ciam_protect, :floodmortality, 0.01)

# Wetland params
setparameter(m2, :ciam_protect, :wbvm, 0.376)
setparameter(m2, :ciam_protect, :wetlandarea, 30.)

# Slr params
setparameter(m2, :ciam_protect, :lslr, lslr)
setparameter(m2, :ciam_protect, :slr10, 0.1)
setparameter(m2, :ciam_protect, :slr100, 0.2)
setparameter(m2, :ciam_protect, :slr1000, 0.3)
setparameter(m2, :ciam_protect, :slr10000, 3.3)


# Sanity check for model outcomes
m2[:ciam_protect, :H10]
m2[:ciam_protect, :H100]
m2[:ciam_protect, :H1000]
m2[:ciam_protect, :H10000]

m2[:ciam_protect, :ProtectCost]
m2[:ciam_protect, :ProtectLevel]
m2[:ciam_protect, :AdaptationOption]
m2[:ciam_protect, :AdaptationLevel]
m2[:ciam_protect, :AdaptationCost]


# Evaluate Model Results/Outcomes
runtimes = ones(1.:100.)
runtimes2 = ones(1.:100.)
for i in collect(1:100)
    runtimes[i] = (@timed run(m))[2]
    runtimes2[i] = (@timed run(m2))[2]
end

mean(runtimes)
mean(runtimes2)