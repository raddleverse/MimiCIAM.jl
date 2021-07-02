using MimiCIAM
using CSV

"""
    write_init_file(run_name::String, outputdir::String, init_settings::Dict)

Write the init.csv file for a specificied `run_name` into `outputdir` using init_settings
found in `init_settings`.  Note the file will be named `<run_name>_init.csv`
"""
function write_init_file(run_name::String, outputdir::String, init_settings::Dict)
    textheader="run_name,lslr,subset,ssp,ssp_simplified\n"
    textstr = "$(run_name),$(init_settings[:lslrfile]),$(init_settings[:subset]),$(init_settings[:ssp]),$(init_settings[:ssp_simplified])"
    txtfile = open(joinpath(outputdir, init_settings[:init_filename]),"w") do io
        write(io,textheader)
        write(io,textstr)
    end
end

"""
    write_output_files(m, outputdir::String, run_name::String)

Write three output files for run `run_name` of model `m` into` outputdir`. These
files include:
- Writing out ciam `subsegs = seg` file for `run_name` to directory `outputdir`
- Writing out ciam `subsegs = global` file for run `run_name` to directory `outputdir`
- Writing out optimal costs file for run `run_name` to directory `outputdir`
"""
function write_output_files(m, outputdir::String, run_name::String)

    # write out the results
    println("Writing out ciam `subsegs = seg` file for run $(run_name) to directory $(outputdir)")
    MimiCIAM.write_ciam(m; outputdir = outputdir, runname = run_name, sumsegs="seg", varnames=false)
    println("Writing out ciam `subsegs = global` file for run $(run_name) to directory $(outputdir)")
    MimiCIAM.write_ciam(m; outputdir = outputdir, runname = run_name, sumsegs="global", varnames=false)
    println("Writing out optimal costs file for run $(run_name) to directory $(outputdir)")
    MimiCIAM.write_optimal_costs(m; outputdir = outputdir, runname = run_name)

end
