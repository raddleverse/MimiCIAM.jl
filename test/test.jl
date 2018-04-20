# Tests For CIAM Model Comparison
# TODO update comments
# TODO fix result overwrite issues
#------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------
using DataFrames
using Plots
using StatPlots
gr()

function compare_outputs(A, B, header, metadata, file, tol=1e-5)

    # Compute pct diff
    if length(A)==length(B)
        pctdiff = abs.((B - A)./B) 
    else
        return nothing
    end
   
    # Parse metadata into string - scenario, segment, category
    meta = string(metadata[1])
    for i in 2:length(metadata)
        meta = string(meta, ",", metadata[i])
    end

    
    # Rank according to priority on 1-5 scale (1 worst 5 best) and write to file
    open(file,"w") do f

        println(f,header)    

        for i in 1:length(pctdiff)
            if pctdiff[i] <= tol
                rank = 6
            elseif pctdiff[i] > tol && pctdiff[i] <= .01
                rank = 5
            elseif pctdiff[i] > .01 && pctdiff[i] <= .1
                rank = 4
            elseif pctdiff[i] > .1 && pctdiff[i] <= .5
                rank = 3
            elseif pctdiff[i] > .5 && pctdiff[i] <= 1
                rank = 2
            else
                rank = 1
            end

            println(f,"$(meta),$(i),$(rank),$(pctdiff[i]*100)")
        end

    end

end

# Function to run CIAM for given parameters and segment-country mapping
# Some params are currently hard-coded
function run_model(params, xsc)
    m = Model()
    setindex(m, :time, 20)
    setindex(m, :adaptPers, 5)  # Todo figure out way not to hardcode these
    setindex(m, :regions, xsc[2])
    setindex(m, :segments, xsc[3])

    addcomponent(m, ciam)
    setparameter(m, :ciam, :xsc, xsc[1])
    setleftoverparameters(m, params)

    run(m)

    return m
end

# Wrapper for importing model data. Not generalizable; hard-coded names; WIP
# datadir - data directory
# lslfile - filename for lsl (string)
# xscfile - filename for segment-country mapping (string)
function import_model_data(datadir, lslfile, xscfile, subset)
    # Process main and lsl params
    params = load_ciam_params(datadir)
     
    # Process XSC
    xsc = prepxsc(datadir, xscfile, subset)

    parse_ciam_params!(params, xsc[2], xsc[3])
    preplsl!(datadir, lslfile, subset, params)

    return(params, xsc)

end

# Function to import comparison data and store in DataFrame
# WIP - to do -- unique header or detect header automatically
function import_comparison_data(datadir, file)
    data = readtable(abspath(datadir,file))
    data = DataFrame(data)
    return(data)
end

function make_plots(plotlist,title,resultsdir)
    if length(plotlist)==3
         p = plot(plotlist[1],plotlist[2],plotlist[3],layout=(3,1),legend=(:bottom))
    elseif length(plotlist)==4
        p = plot(plotlist[1],plotlist[2],plotlist[3],plotlist[4], layout=(2,2),legend=(:bottom))
    elseif length(plotlist)==5
        p = plot(plotlist[1],plotlist[2],plotlist[3],plotlist[4],plotlist[5],plot(0), layout=(2,3),legend=(:bottom))
    else
        p = plot(plotlist[1], legend=(:bottom))
    end
    savefig(p, joinpath(resultsdir,"$(title)_plot.pdf"))
end

function line_plot(gams, jl, title)

    p = plot(gams, label="GAMS",yaxis=("billion2010USD"))
    plot!(p, jl, title = title, linestyle = :dash, label = "Julia")

    return p

end

# Function to write out model results to CSV file
# m - an instance of the model
# RCP - string for RCP we're using; todo make automatic from lsl file
# xsc - segment-region dictionaries
# outputdir, outfile - where to write results to, relative to test folder
# sumsegs - whether to sum across segments
function write_results(m, rcp, outputdir, xsc, sumsegs = true)

    meta_output = load_meta()
    outfile = "results_$(rcp).csv"
    header = meta_output[1]
    vardict = meta_output[2]
    protectdict = meta_output[3]
    retreatdict = meta_output[4]

    segmap = xsc[6]

    vars = [k for k in keys(vardict)]

    open(joinpath(outputdir,outfile),"w") do g 
        write(g, "$(header)\n")
    end

    for v in vars # Iterate variables
        d = m[:ciam, parse(v)]
        
        level = vardict[v][1]
        costtype = vardict[v][2]

        s_arr = [segmap[i] for i in 1:size(d,1)]    # List of segments
        
        if sumsegs
            func = x -> sum(x)
            n = 1
            s_arr = ["sum$(size(d,1))segs"]
        else
            func = x -> identity(x)
            n = size(d,1)
        end

        for i in 1:size(d,2) # Iterate times (t = 1, 2, ...)
            t = repeat([i], outer=n)
            rcp_arr = repeat([rcp], outer= n)
            cost_arr = repeat([costtype], outer = n)

            if typeof(d)==Array{Float64,2}
                val = d[:,i] / 10 * .001
                v_arr = func(val)

                if level=="protect"
                    for j in 1:4
                        p = protectdict[j]
                        lev_arr = repeat([p], outer = n)

                        outarr= [rcp_arr lev_arr s_arr cost_arr t v_arr]
                            # Write results to csv 
             
                        open(joinpath(outputdir,outfile),"a") do g 
                            writecsv(g, outarr)
                        end

                    end
                elseif level=="retreat"
                    for j in 1:5
                        q = retreatdict[j]
                        lev_arr = repeat([q], outer = n)

                        outarr= [rcp_arr lev_arr s_arr cost_arr t v_arr]  
            
                        open(joinpath(outputdir,outfile),"a") do g 
                            writecsv(g, outarr)
                        end
                    end

                else
                    lev_arr = repeat([level], outer = n)
                    outarr= [rcp_arr lev_arr s_arr cost_arr t v_arr]

                    open(joinpath(outputdir,outfile),"a") do g 
                        writecsv(g, outarr)
                    end

                end

            elseif typeof(d)==Array{Float64,3}
                for j in 1:size(d,3)
                    val = d[:,i,j] / 10 * .001
                    v_arr = func(val)

                    if level=="protect"
                        p = protectdict[j]
                        lev_arr = repeat([p], outer = n)

                    elseif level=="retreat"
                        q = retreatdict[j]
                        lev_arr = repeat([q], outer = n)
                    else
                        lev_arr = repeat([level], outer = n)
                    end

                    outarr= [rcp_arr lev_arr s_arr cost_arr t v_arr]
                 
                    open(joinpath(outputdir,outfile),"a") do g 
                        writecsv(g, outarr)
                    end

                end
            end

        end
    end

end

function load_meta()
    metadir = "../data/meta"

    # Read header and mappings
    header = readstring(open(joinpath(metadir,"header.txt")))

    varnames = readlines(open(joinpath(metadir,"variablenames.csv")))
    vardict = Dict{Any,Any}(split(varnames[i],',')[1] => (split(varnames[i],',')[2],split(varnames[i],',')[3]) for i in 1:length(varnames))

    protect = readlines(open(joinpath(metadir, "protectlevels.csv")))
    protectdict = Dict{Any,Any}(parse(Int,split(protect[i],',')[1]) => split(protect[i],',')[2] for i in 1:length(protect))

    retreat = readlines(open(joinpath(metadir, "retreatlevels.csv")))
    retreatdict = Dict{Any,Any}(parse(Int,split(retreat[i],',')[1]) => split(retreat[i],',')[2] for i in 1:length(retreat))

    return(header,vardict, protectdict, retreatdict)
end


# Function to run tests for model and make plots
# WIP - hardcoded strings 
# datadir - directory where data is located; 
# resultsdir - output directory
# gamsdata,jldata - location of comparison data from GAMS/Julia (relative to datadir)
# rcp - string rcp value to test
function run_tests(datadir, gamsfile, resultsdir, lslfile, subset, rcp, model=true, conv_factor = 1e-4)
    # Import model data
    modelparams = import_model_data(datadir, lslfile,"xsc.csv", subset)
    params = modelparams[1]
    xsc = modelparams[2]

    # Import GAMS Data
    gamsdata = import_comparison_data(datadir, gamsfile)

    # Load Metadata
    metavars = load_meta()
    
    # Run model if specified
    if model

        m = run_model(params, xsc)

        if (length(subset)>5 || subset==false)
            sum = true
        else
            sum = false
        end

        write_results(m, rcp, resultsdir, xsc, sum)
    end

    # Compare and ouptut results
    jldata = import_comparison_data(resultsdir, "results_$(rcp).csv") # TODO - distinct outputs for each model run
    levels = readlines(open("../data/meta/levels.csv"))                  
    variables = unique([v[2] for v in values(metavars[2])])
    segments = unique(jldata[:seg])

    for j in 1:length(levels)
        plotlist = []
        for k in 1:length(variables)

            # Skip incompatible combinations
            if (contains(levels[j], "retreat")|| levels[j]=="noAdaptation") && variables[k]=="protection"
                continue
            elseif contains(levels[j], "protect") && (variables[k]=="inundation" || variables[k]=="relocation")
                continue
            elseif (levels[j]=="optimalFixed") && (variables[k] != "total") # TODO update this 
                continue
            else
        
                if sum
                    A = jldata[ (jldata[:level] .== levels[j]) .& (jldata[:costtype] .== variables[k]), :value]
                    B = gamsdata[ (gamsdata[:level] .== levels[j]) .& (gamsdata[:costtype] .== variables[k]), :value]

                     # Make plots
                     title = string(unique(jldata[:seg])[1], variables[k])
                     p = line_plot(B, A, title)
                     push!(plotlist, p)
                else
                    for s in segments
                        A = jldata[ (jldata[:level] .== levels[j]) .& (jldata[:costtype] .== variables[k]) .& (jldata[:seg] .== s), :value] 
                        B = gamsdata[ (gamsdata[:level] .== levels[j]) .& (gamsdata[:costtype] .== variables[k]) .& (gamsdata[:seg].==s), :value]

                        # Make plots
                        title = string(s,variables[k])
                        p = line_plot(B, A, title)
                        push!(plotlist, p)
                    end

                end        
            end
                
        end
        make_plots(plotlist,"$(rcp)_$(levels[j])", resultsdir)
    end

    if model
          return m 
    else 
        return nothing
    end

end


