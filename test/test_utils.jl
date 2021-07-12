using CSV

##==============================================================================
## Write Comparison Files

function write_MimiCIAM_comparison_files(outputdir; subset = false)

    isdir(outputdir) || mkpath(outputdir)

    ##==============================================================================
    ## ctrl+noConstrFix: This case is run with a modified slrcost component held in
    ## slrcost_GAMSmatch.jl, which is taken care of in the `get_model` step with the
    ## GAMS match arg and removes the block that disallows height reductions

    run_name = "ctrl+noConstrFix"

    init_settings = Dict(
        :init_filename   => string("$run_name", "_init.csv"),
        :lslrfile        => "lsl_rcp85_p50.csv",
        :subset          => subset,
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
    MimiCIAM.write_init_file(run_name, outputdir, init_settings)

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

    MimiCIAM.write_output_files(m, outputdir, run_name)

    ##==============================================================================
    ##  Ctrl Case

    run_name = "ctrl"

    init_settings = Dict(
        :init_filename   => string("$run_name", "_init.csv"),
        :lslrfile        => "lsl_rcp85_p50.csv",
        :subset          => subset,
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

    MimiCIAM.write_init_file(run_name, outputdir, init_settings)

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

    MimiCIAM.write_output_files(m, outputdir, run_name)

    ##==============================================================================
    ##  baseline+updated GDP/POP via SSP5. but can be any of 1-5

    run_name = "ctrl+SSP5"

    init_settings = Dict(
        :init_filename   => string("$run_name", "_init.csv"),
        :lslrfile        => "lsl_rcp85_p50.csv",
        :subset          => subset,
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

    MimiCIAM.write_init_file(run_name, outputdir, init_settings)

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

    MimiCIAM.write_output_files(m, outputdir, run_name)

end
