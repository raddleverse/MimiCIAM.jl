
##==============================================================================
## baseline_comparisons.jl
## baseline comparisons against the GAMS results
##
## Tony Wong (aewsma@rit.edu)
##==============================================================================

##==============================================================================
## baseline cases
noRetreat = false
allowMaintain = false
fixed = true
subset = false
#subset = "segments_test.csv"
#subset = "sub1names.csv"

## no SLR case and RCP8.5 SLR cases (median and 5th/95th percentiles)
#for b in ["lsl_rcp85_p50.csv"]#"lsl_rcp0_p50.csv","lsl_rcp85_p50.csv","lsl_rcp85_p5.csv","lsl_rcp85_p95.csv"]
pop = 0
SSP = 0
SSP_simp = 2 # won't matter, just to get defaults for the popdens_seg_jones and _merkens arrays.
runname="ctrl+noConstrFix"
b= "lsl_rcp85_p50.csv"
textheader="run_name,lslr,subset,ssp,ssp_simplified\n"
textstr = "base,$(b),$(subset),$(SSP[1]),$(SSP_simp)"
txtfile=open("../data/batch/init.txt","w") do io
    write(io,textheader)
    write(io,textstr)
end
m = MimiCIAM.get_model(t=15,initfile="../data/batch/init.txt",fixed=fixed,noRetreat=noRetreat,allowMaintain=allowMaintain)
run(m)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="global", varnames=false)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="seg", varnames=false)
MimiCIAM.write_optimal_costs(m; runname=runname)
    # need the breakdown of optimal costs since they're different for each segment
    #  (different adaptation strategies)
#    if b in ["lsl_rcp0_p50.csv","lsl_rcp85_p50.csv"]
        #MimiCIAM.write_optimal_costs(m; runname=runname)
#    end
#end
##==============================================================================


##==============================================================================
## changing things

##=====================================
## here, need to uncomment the Hprev > H fix in slrcost.jl, re-start Julia and re-run.
## this should be uncommented and used from here forward, which is why it is not an optional parameter
using Mimi, MimiCIAM, Query, RData, StatsBase, CSV, DataFrames, NetCDF
cd("/Users/aewsma/codes/CIAM_adaptation_regimes/ciam-code/src")
noRetreat = false
allowMaintain = false
fixed = true
subset = false
pop = 0
SSP = 0
SSP_simp = 2 # won't matter, just to get defaults for the popdens_seg_jones and _merkens arrays.
runname="ctrl"
b= "lsl_rcp85_p50.csv"
textheader="run_name,lslr,subset,ssp,ssp_simplified\n"
textstr = "base,$(b),$(subset),$(SSP),$(SSP_simp)"
txtfile=open("../data/batch/init.txt","w") do io
    write(io,textheader)
    write(io,textstr)
end
m = MimiCIAM.get_model(t=20,initfile="../data/batch/init.txt",fixed=fixed,noRetreat=noRetreat,allowMaintain=allowMaintain)
run(m)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="global", varnames=false)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="seg", varnames=false)
MimiCIAM.write_optimal_costs(m; runname=runname)
##=====================================

##=====================================
## at this point, with updated SLR, need to modify `at.csv` to exclude the adaptation
## period starting at t=19 because that is outside the time horizon of the Vega-Westhoff et al 2020

## baseline+updated GDP/POP via SSP5. but can be any of 1-5
SSP = "IIASAGDP_SSP5_v9_130219"
SSP_simp = "5"
runname="ctrl+SSP5"
b= "lsl_rcp85_p50.csv"
textheader="run_name,lslr,subset,ssp,ssp_simplified\n"
textstr = "base,$(b),$(subset),$(SSP),$(SSP_simp)"
txtfile=open("../data/batch/init.txt","w") do io
    write(io,textheader)
    write(io,textstr)
end
m = MimiCIAM.get_model(t=15,initfile="../data/batch/init.txt",fixed=fixed,noRetreat=noRetreat,allowMaintain=allowMaintain)
run(m)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="global", varnames=false)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="seg", varnames=false)
MimiCIAM.write_optimal_costs(m; runname=runname)


## baseline+updated population density from Jones and O'Neill (2016)
pop = 1 # 0 = default old stuff, 1 = Jones and O'Neill (2016), 2 = Merkens et al (2016)
SSP = "IIASAGDP_SSP5_v9_130219"
SSP_simp = "5" # based on SSPs, so need to use an SSP
runname="ctrl+SSP5+popJones"
b= "lsl_rcp85_p50.csv"
textheader="run_name,lslr,subset,ssp,ssp_simplified\n"
textstr = "base,$(b),$(subset),$(SSP),$(SSP_simp)"
txtfile=open("../data/batch/init.txt","w") do io
    write(io,textheader)
    write(io,textstr)
end
m = MimiCIAM.get_model(t=15,initfile="../data/batch/init.txt",fixed=fixed,noRetreat=noRetreat,allowMaintain=allowMaintain, popinput=pop)
update_param!(m, :ssp, parse(Int32, SSP_simp))
run(m)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="global", varnames=false)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="seg", varnames=false)
MimiCIAM.write_optimal_costs(m; runname=runname)


## baseline+updated SLR RCP8.5
include("brickLSL.jl")
rcp = 85
brickfile = "../data/lslr/BRICK_projections.RData"
runname="ctrl+BRICKLSL85"
pop = 0
SSP = 0  # to reset SSP
SSP_simp = 2 # won't matter, just to get defaults for the popdens_seg_jones and _merkens arrays.
textheader="run_name,lslr,subset,ssp,ssp_simplified\n"
textstr = "base,$(b),$(subset),$(SSP),$(SSP_simp)"
txtfile=open("../data/batch/init.txt","w") do io
    write(io,textheader)
    write(io,textstr)
end
lsl = brick_lsl(rcp,segIDs,brickfile,1,50,50,2010,2150,10,false) # end_year here only for picking median SLR ensemble member
lslr=lsl[1]
gmsl=lsl[2]
ensInds=lsl[3] # Indices of original BRICK array
m = MimiCIAM.get_model(t=15,initfile="../data/batch/init.txt",fixed=fixed,noRetreat=noRetreat,allowMaintain=allowMaintain)
update_param!(m,:lslr,lslr[1,:,:])
run(m)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="global", varnames=false)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="seg", varnames=false)
MimiCIAM.write_optimal_costs(m; runname=runname)

## baseline+updated SLR RCP2.6
rcp = 26
brickfile = "../data/lslr/BRICK_projections.RData"
runname="ctrl+BRICKLSL26"
textheader="run_name,lslr,subset,ssp,ssp_simplified\n"
textstr = "base,$(b),$(subset),$(SSP),$(SSP_simp)"
txtfile=open("../data/batch/init.txt","w") do io
    write(io,textheader)
    write(io,textstr)
end
lsl = brick_lsl(rcp,segIDs,brickfile,1,50,50,2010,2150,10,false) # end_year here only for picking median SLR ensemble member
lslr=lsl[1]
gmsl=lsl[2]
ensInds=lsl[3] # Indices of original BRICK array
m = MimiCIAM.get_model(t=15,initfile="../data/batch/init.txt",fixed=fixed,noRetreat=noRetreat,allowMaintain=allowMaintain)
update_param!(m,:lslr,lslr[1,:,:])
run(m)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="global", varnames=false)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="seg", varnames=false)
MimiCIAM.write_optimal_costs(m; runname=runname)


## baseline+updated SLR (RCP8.5)+updated GPD/POP via SSP5+updated population density from Jones and O'Neill (2016)
rcp = 85
brickfile = "../data/lslr/BRICK_projections.RData"
runname="ctrl+SSP5+popJones+BRICKLSL85"
pop = 1
SSP = "IIASAGDP_SSP5_v9_130219"
SSP_simp = "5" #SSP = "IIASAGDP_SSP5_v9_130219"
lsl = brick_lsl(rcp,segIDs,brickfile,1,50,50,2010,2150,10,false) # end_year here only for picking median SLR ensemble member
lslr=lsl[1]
gmsl=lsl[2]
ensInds=lsl[3] # Indices of original BRICK array
m = MimiCIAM.get_model(t=15,initfile="../data/batch/init.txt",fixed=fixed,noRetreat=noRetreat,allowMaintain=allowMaintain, popinput=pop)
update_param!(m, :lslr, lslr[1,:,:])
update_param!(m, :ssp, parse(Int32, SSP_simp))
run(m)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="global", varnames=false)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="seg", varnames=false)
MimiCIAM.write_optimal_costs(m; runname=runname)

## same, but ensemble member giving the 5th percentile of GMSL used
rcp = 85
brickfile = "../data/lslr/BRICK_projections.RData"
runname="ctrl+SSP5+popJones+BRICKLSL85_p05"
pop = 1
SSP = "IIASAGDP_SSP5_v9_130219"
SSP_simp = "5" #SSP = "IIASAGDP_SSP5_v9_130219"
lsl = brick_lsl(rcp,segIDs,brickfile,1,5,5,2010,2150,10,false) # end_year here only for picking median SLR ensemble member
lslr=lsl[1]
gmsl=lsl[2]
ensInds=lsl[3] # Indices of original BRICK array
m = MimiCIAM.get_model(t=15,initfile="../data/batch/init.txt",fixed=fixed,noRetreat=noRetreat,allowMaintain=allowMaintain, popinput=pop)
update_param!(m, :lslr, lslr[1,:,:])
update_param!(m, :ssp, parse(Int32, SSP_simp))
run(m)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="global", varnames=false)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="seg", varnames=false)
MimiCIAM.write_optimal_costs(m; runname=runname)

rcp = 85
brickfile = "../data/lslr/BRICK_projections.RData"
runname="ctrl+SSP5+popJones+BRICKLSL85_p95"
pop = 1
SSP = "IIASAGDP_SSP5_v9_130219"
SSP_simp = "5" #SSP = "IIASAGDP_SSP5_v9_130219"
lsl = brick_lsl(rcp,segIDs,brickfile,1,95,95,2010,2150,10,false) # end_year here only for picking median SLR ensemble member
lslr=lsl[1]
gmsl=lsl[2]
ensInds=lsl[3] # Indices of original BRICK array
m = MimiCIAM.get_model(t=15,initfile="../data/batch/init.txt",fixed=fixed,noRetreat=noRetreat,allowMaintain=allowMaintain, popinput=pop)
update_param!(m, :lslr, lslr[1,:,:])
update_param!(m, :ssp, parse(Int32, SSP_simp))
run(m)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="global", varnames=false)
MimiCIAM.write_ciam(m; runname=runname, sumsegs="seg", varnames=false)
MimiCIAM.write_optimal_costs(m; runname=runname)


##==============================================================================


##==============================================================================
## End
##==============================================================================