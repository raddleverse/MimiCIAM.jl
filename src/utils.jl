using MimiCIAM
using CSV

## Writing helper funcitons

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


# Helper functions for slrcost component
function growthrate(x1, x2)
    epsilon = 1.e-9
    return (x2 / (x1 + epsilon) - 1.)
end

function calcCoastArea(areaparams, var)

    area = (areaparams[1]*max(0,min(0.5,var-0))
    +(areaparams[1]+areaparams[2])/2*max(0,min(1,var-0.5))
    +areaparams[2]*max(0,min(0.5,var-1.5))
    +areaparams[3]*max(0,min(1,var-2))
    +areaparams[4]*max(0,min(1,var-3))
    +areaparams[5]*max(0,min(1,var-4))
    +areaparams[6]*max(0,min(1,var-5))
    +areaparams[7]*max(0,min(1,var-6))
    +areaparams[8]*max(0,min(1,var-7))
    +areaparams[9]*max(0,min(1,var-8))
    +areaparams[10]*max(0,min(1,var-9))
    +areaparams[11]*max(0,min(1,var-10))
    +areaparams[12]*max(0,min(1,var-11))
    +areaparams[13]*max(0,min(1,var-12))
    +areaparams[14]*max(0,min(1,var-13))
    +areaparams[15]*max(0,var-14))

    return area

end

function localrate(lslr1, lslr2, tstep)
    return max(0, (lslr2 - lslr1))/tstep
end

function getsegments(rgn_name, xsc)

    segs = collect(keys(filter( (k,v) -> v[1]==rgn_name,xsc)))
    return segs

end

function getregion(seg_ind, xsc)
    rgn = xsc[seg_ind][1]

    return rgn
end

function isgreenland(seg_ind, xsc)
    greenland = xsc[seg_ind][2]
    return greenland
end

function isisland(seg_ind, xsc)
    island = xsc[seg_ind][3]
    return island
end

function calcHorR(option, level, lslrPlan, surgeExpLevels, adaptOptions)

    ind = findind(level, adaptOptions)

    if option==-1 && level ==10
        # Protect height differs from retreat radius only in case of 10 yr surge exposure
        H = max(0, lslrPlan + surgeExpLevels[ind] / 2)
        return H
    elseif level==0 # Maintain existing defenses
        H_R=0
        return H_R
    else
        H_R = max(0, lslrPlan + surgeExpLevels[ind])
        return H_R
    end
end

function pos(x)
    return max(0,x)
end