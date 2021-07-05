using MimiCIAM

## Setup

outputdir = joinpath(@__DIR__, "..", "output", "main")
isdir(outputdir) || mkpath(outputdir)

## A Basic Run

m = MimiCIAM.get_model()
run(m)
MimiCIAM.write_output_files(m, outputdir, "default") # writes three files

## Write out optimal costs for all SLR options

lsl_files = [
    "BRICKsneasy-lsl_rcp85_p5.csv",
    "BRICKsneasy-lsl_rcp85_p10.csv",
    "BRICKsneasy-lsl_rcp85_p50.csv",
    "BRICKsneasy-lsl_rcp85_p90.csv",
    "BRICKsneasy-lsl_rcp85_p95.csv",
    "lsl_rcp0_p50.csv",
    "lsl_rcp45_p50.csv",
    "lsl_rcp85_p5.csv",
    "lsl_rcp85_p50.csv",
    "lsl_rcp85_p95.csv"
]

for file in lsl_files

    run_name = string(split(file, ".")[1])
    run_outputdir = joinpath(outputdir, run_name)
    isdir(run_outputdir) || mkpath(run_outputdir)

    init_settings = Dict(
        :init_filename   => string("$run_name", "_init.csv"),
        :lslrfile        => file,
        :subset          => false,
        :ssp             => "IIASAGDP_SSP5_v9_130219",
        :ssp_simplified  => 5
    )

    MimiCIAM.write_init_file(run_name, run_outputdir, init_settings)
    m = MimiCIAM.get_model(initfile = joinpath(run_outputdir, init_settings[:init_filename]))
    run(m)

    MimiCIAM.write_optimal_costs(m; outputdir = run_outputdir, runname = run_name)
    # MimiCIAM.write_output_files(m; outputdir = outputdir, run_name = run_name) # writes three files

end
