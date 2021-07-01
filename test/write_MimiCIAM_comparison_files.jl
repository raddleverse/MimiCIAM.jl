using CSV
using MimiCIAM
using Query
using RData
using StatsBase
using CSV
using DataFrames
using NetCDF

##==============================================================================
## Setup and Helper Functions

outputdir = joinpath(@__DIR__, "..", "output", "results-jl")
isdir(outputdir) || mkdir(outputdir)

function write_init_file(run_name::String, outputdir::String, init_settings::Dict)
    textheader="run_name,lslr,subset,ssp,ssp_simplified\n"
    textstr = "$(run_name),$(init_settings[:lslrfile]),$(init_settings[:subset]),$(init_settings[:ssp]),$(init_settings[:ssp_simplified])"
    txtfile = open(joinpath(outputdir, init_settings[:init_filename]),"w") do io
        write(io,textheader)
        write(io,textstr)
    end
end

function write_MimiCIAM_comparison_files(m, outputdir::String, run_name::String)

    # write out the results
    println("Writing out ciam `subsegs = seg` file for run $(run_name) to directory $(outputdir)")
    MimiCIAM.write_ciam(m; outputdir = outputdir, runname = run_name, sumsegs="seg", varnames=false)
    println("Writing out ciam `subsegs = global` file for run $(run_name) to directory $(outputdir)")
    MimiCIAM.write_ciam(m; outputdir = outputdir, runname = run_name, sumsegs="global", varnames=false)
    println("Writing out optimal costs file for run $(run_name) to directory $(outputdir)")
    MimiCIAM.write_optimal_costs(m; outputdir = outputdir, runname = run_name)

end

##==============================================================================
## Write Comparison Files

##==============================================================================
## ctrl+noConstrFix: This case is run with a modified slrcost component held in 
## slrcost_GAMSmatch.jl, which is taken care of in the `get_model` step with the 
## GAMS match arg and removes the block that disallows height reductions

run_name = "ctrl+noConstrFix"

init_settings = Dict(
    :init_filename   => string("$run_name", "_init.csv"),
    :lslrfile        => "lsl_rcp85_p50.csv",
    :subset          => false,
    :ssp             => 0,
    :ssp_simplified  => 2 # won't matter, just to get defaults for the popdens_seg_jones and _merkens arrays.
)

model_settings = Dict(
    :fixed          => true,
    :t              => 20,
    :noRetreat      => false,
    :allowMaintain  => false,
    :popinput       => 0,
    :GAMSmatch      => true
)

# write files
write_init_file(run_name, outputdir, init_settings)

m = MimiCIAM.get_model(
    initfile        = joinpath(outputdir, init_settings[:init_filename]),
    fixed           = model_settings[:fixed],
    t               = model_settings[:t], 
    noRetreat       = model_settings[:noRetreat],
    allowMaintain   = model_settings[:allowMaintain],
    popinput        = model_settings[:popinput],
    GAMSmatch       = model_settings[:GAMSmatch]
)
run(m)

write_MimiCIAM_comparison_files(run_name, outputdir, run_name)

##==============================================================================
##  Ctrl Case

run_name = "ctrl"

init_settings = Dict(
    :init_filename   => string("$run_name", "_init.csv"),
    :lslrfile        => "lsl_rcp85_p50.csv",
    :subset          => false,
    :ssp             => 0,
    :ssp_simplified  => 2 # won't matter, just to get defaults for the popdens_seg_jones and _merkens arrays.
)

model_settings = Dict(
    :fixed          => true,
    :t              => 15,
    :noRetreat      => false,
    :allowMaintain  => false,
    :popinput       => 0,
    :GAMSmatch      => false
)

write_init_file(run_name, outputdir, init_settings)

m = MimiCIAM.get_model(
    initfile        = joinpath(outputdir, init_settings[:init_filename]),
    fixed           = model_settings[:fixed],
    t               = model_settings[:t], 
    noRetreat       = model_settings[:noRetreat],
    allowMaintain   = model_settings[:allowMaintain],
    popinput        = model_settings[:popinput],
    GAMSmatch       = model_settings[:GAMSmatch]
)
run(m)

write_MimiCIAM_comparison_files(m, outputdir, run_name)

##==============================================================================
##  baseline+updated GDP/POP via SSP5. but can be any of 1-5

run_name = "ctrl+SSP5"

init_settings = Dict(
    :init_filename   => string("$run_name", "_init.csv"),
    :lslrfile        => "lsl_rcp85_p50.csv",
    :subset          => false,
    :ssp             => "IIASAGDP_SSP5_v9_130219",
    :ssp_simplified  => 5
)

model_settings = Dict(
    :fixed          => true,
    :t              => 15,
    :noRetreat      => false,
    :allowMaintain  => false,
    :popinput       => 0,
    :GAMSmatch      => false
)

write_init_file(run_name, outputdir, init_settings)

m = MimiCIAM.get_model(
    initfile        = joinpath(outputdir, init_settings[:init_filename]),
    fixed           = model_settings[:fixed],
    t               = model_settings[:t], 
    noRetreat       = model_settings[:noRetreat],
    allowMaintain   = model_settings[:allowMaintain],
    popinput        = model_settings[:popinput],
    GAMSmatch       = model_settings[:GAMSmatch]
)
run(m)

write_MimiCIAM_comparison_files(m, outputdir, run_name)

##==============================================================================
##  baseline+updated population density from Jones and O'Neill (2016)

run_name = "ctrl+SSP5+popJones"

init_settings = Dict(
    :init_filename   => string("$run_name", "_init.csv"),
    :lslrfile        => "lsl_rcp85_p50.csv",
    :subset          => false,
    :ssp             => "IIASAGDP_SSP5_v9_130219",
    :ssp_simplified  => "5" # based on SSPs, so need to use an SSP
)

model_settings = Dict(
    :fixed          => true,
    :t              => 15,
    :noRetreat      => false,
    :allowMaintain  => false,
    :popinput       => 1, # 0 = default old stuff, 1 = Jones and O'Neill (2016), 2 = Merkens et al (2016)
    :GAMSmatch      => false
)

write_init_file(run_name, outputdir, init_settings)

m = MimiCIAM.get_model(
    initfile        = joinpath(outputdir, init_settings[:init_filename]),
    fixed           = model_settings[:fixed],
    t               = model_settings[:t], 
    noRetreat       = model_settings[:noRetreat],
    allowMaintain   = model_settings[:allowMaintain],
    popinput        = model_settings[:popinput],
    GAMSmatch       = model_settings[:GAMSmatch]
)
update_param!(m, :ssp, parse(Int32, SSP_simp))
run(m)

write_MimiCIAM_comparison_files(m, outputdir, run_name)
