using Mimi
using MimiCIAM
using CSV
using RData
using DataFrames 

function write_init_file(run_name::String, outputdir::String, init_settings::Dict)
    textheader="run_name,lslr,subset,ssp,ssp_simplified\n"
    textstr = "$(run_name),$(init_settings[:lslrfile]),$(init_settings[:subset]),$(init_settings[:ssp]),$(init_settings[:ssp_simplified])"
    txtfile = open(joinpath(outputdir, init_settings[:init_filename]),"w") do io
        write(io,textheader)
        write(io,textstr)
    end
end

function write_output_files(m, outputdir::String, run_name::String)
    println("Writing out ciam `subsegs = seg` file for run $(run_name) to directory $(outputdir)")
    MimiCIAM.write_ciam(m; outputdir = outputdir, runname = run_name, sumsegs="seg", varnames=false)
    println("Writing out ciam `subsegs = global` file for run $(run_name) to directory $(outputdir)")
    MimiCIAM.write_ciam(m; outputdir = outputdir, runname = run_name, sumsegs="global", varnames=false)
    println("Writing out optimal costs file for run $(run_name) to directory $(outputdir)")
    MimiCIAM.write_optimal_costs(m; outputdir = outputdir, runname = run_name)
end
##==============================================================================
## Setup

# outputdir = joinpath(@__DIR__, "..", "output", "baseline_comparisons")
outputdir = "/Users/lisarennels/JuliaProjects/CIAMPaper/local-data/model-outputs-preLisa"
isdir(outputdir) || mkpath(outputdir)
##==============================================================================
## Run Comparisons

##==============================================================================
## ctrl+noConstrFix

## To run this one the `at` file should go to 19 and the height section should
## be commented out of the slrcost component

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
    :popinput       => 0
    # :GAMSmatch      => true
)

write_init_file(run_name, outputdir, init_settings)

m = MimiCIAM.get_model(
    initfile        = joinpath(outputdir, init_settings[:init_filename]),
    fixed           = model_settings[:fixed],
    t               = model_settings[:t], 
    noRetreat       = model_settings[:noRetreat],
    allowMaintain   = model_settings[:allowMaintain],
    popinput        = model_settings[:popinput]
    # GAMSmatch       = model_settings[:GAMSmatch]
)
run(m)

write_output_files(m, outputdir, run_name)

##==============================================================================
##  Ctrl Case

## To run this one the `at` file should go to 15 and the height section should NOT
## be commented out of the slrcost component

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
    :popinput       => 0
    # :GAMSmatch      => false
)

write_init_file(run_name, outputdir, init_settings)

m = MimiCIAM.get_model(
    initfile        = joinpath(outputdir, init_settings[:init_filename]),
    fixed           = model_settings[:fixed],
    t               = model_settings[:t], 
    noRetreat       = model_settings[:noRetreat],
    allowMaintain   = model_settings[:allowMaintain],
    popinput        = model_settings[:popinput]
    # GAMSmatch       = model_settings[:GAMSmatch]
)
run(m)

write_output_files(m, outputdir, run_name)

##==============================================================================
##  baseline+updated GDP/POP via SSP5. but can be any of 1-5

## To run this one the `at` file should go to 15 and the height section should NOT
## be commented out of the slrcost component

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
    :popinput       => 0
    # :GAMSmatch      => false
)

write_init_file(run_name, outputdir, init_settings)

m = MimiCIAM.get_model(
    initfile        = joinpath(outputdir, init_settings[:init_filename]),
    fixed           = model_settings[:fixed],
    t               = model_settings[:t], 
    noRetreat       = model_settings[:noRetreat],
    allowMaintain   = model_settings[:allowMaintain],
    popinput        = model_settings[:popinput]
    # GAMSmatch       = model_settings[:GAMSmatch]
)
run(m)

write_output_files(m, outputdir, run_name)
