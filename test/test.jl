# Testing Suite For CIAM Model Comparison

#------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------
## Function to compare two vector and write results to file
## WiP - still uses hardcoded header string
## A and B are vectors. A is julia and B is GAMS. 
##  metadata - rcp,level,seg,costtype
## file is file to write to
using DataFrames
using Plots
using StatPlots
gr()

function compare_outputs(A, B, metadata, file, tol=1e-5)

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

    if isfile(file)
        option = "a"
        header= false
    else
        option = "w"
        header = true
    end
    
    # Rank according to priority on 1-5 scale (1 worst 5 best) and write to file
    open(file,option) do f

        if header 
            println(f,"rcp,level,seg,costtype,time,rank,pctdiff")    # TODO remove hardcoding
        end
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

## Function to write a summary report from compare_outputs data
## Hard-coded / specific; overwrites previous 
function summary_report(datadir, file, outputdir)
    data = import_comparison_data(datadir, file)

    cases = unique(data[:level])
    costtypes = unique(data[:costtype])
    segments = unique(data[:seg])
    rcps = unique(data[:rcp])

    for i in 1:length(rcps)
        open(joinpath(outputdir, "summary$(rcps[i]).txt"),"w") do f
            for m in 1:length(segments)
                
                println(f,"RCP: $(rcps[i])\tSegment: $(segments[m])")
                
                for j in 1:length(cases)
                    for k in 1:length(costtypes)    
                        df = data[ (data[:rcp].==rcps[i]) .& (data[:seg] .== segments[m]) .& (data[:level] .== cases[j]) .& (d[:costtype].==costtypes[k]) ]
                        if costtypes[k] != "total" && minimum(df[:rank]) <= 3
                            maxdiff = maximum(df[:pctdiff])
                            maxtime = df[ (df[:pctdiff] .== maxdiff) , :time ]
                            println("$(cases[j]): $(costtypes[k]) has deviation >10%; max $(maxdiff)% at t=$(maxtime)")
                        elseif costtypes[k]=="total"
                            mindev = minimum(df[:pctdiff])
                            maxdev = maximum(df[:pctdiff])
                            println("Total: mean rank $(round(mean(df[:rank]))); mindev $(mindev); maxdev $(maxdev)")
                        end

                    end
                end

            end
        end
    end

end

# Function to run CIAM for given parameters and segment-country mapping
# Some params are currently hard-coded
function run_model(params, xsc, model="ciam")
    m = Model()
    setindex(m, :time, 20)
    setindex(m, :adaptPers, 5)  # Todo figure out way not to hardcode these
    setindex(m, :regions, xsc[2])
    setindex(m, :segments, xsc[3])

    addcomponent(m, eval(parse(model)))
    setparameter(m, parse(model), :xsc, xsc[1])
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

    # Parse Main Parameters
    # TODO switch away from using mainparams
    mainparams = Dict{Any,Any}()
    mainparams["data"] = params["data"]
    mainparams["globalparams"] = params["globalparams"]
    mainparams["ypc_usa"] = params["ypc_usa"]
    mainparams["cci"] = params["cci"]
    mainparams["refpopdens"] = params["refpopdens"]
    mainparams["gtapland"] = params["gtapland"]
    mainparams["pop"] = params["pop"]
    mainparams["ypcc"] = params["ypcc"]
    mainparams["at"] = params["at"]
    mainparams["adaptOptions"] = params["adaptoptions"]
    parse_ciam_params!(mainparams, xsc[2], xsc[3])
    preplsl!(datadir, lslfile, subset, mainparams)

    return(mainparams, xsc)

end

# Function to import comparison data and store in DataFrame
# WIP - to do -- unique header or detect header automatically
function import_comparison_data(datadir, file)

    data = readtable(abspath(datadir,file))
    data = DataFrame(data)
    return(data)
end

function make_plots(plotlist,title)
    if length(plotlist)==3
         p = plot(plotlist[1],plotlist[2],plotlist[3],layout=(3,1),legend=(:bottom))
    else
        p = plot(plotlist[1],plotlist[2],plotlist[3],plotlist[4], layout=(2,2),legend=(:bottom))
    end
    savefig(p, "$(title)_plot.pdf")
end

function line_plot(gams, jl, title)
    x = 1:20
    p = plot(x, gams, label="GAMS",yaxis=("billion2010USD"))
    plot!(p, x, jl, title = title, label = "Julia")

    return p

end

# Function to write out model results to CSV file
# m - an instance of the model
# name - model name (e.g. CIAM); string
# meta - model metadata -- e.g. segment name, rcp - from file
# vars - variables to write out; string array
# QUESTION: do we want to a) translate model variable names to results output (e.g. variable name->results dictionary?
#   or b) change results file to use our variable names?
function write_results(m, rcp, outputdir, xsc, truncate = false, sumsegs = true, name = "ciam", outfile = "results.csv")

    meta_output = load_meta(truncate)
    header = meta_output[1]
    vardict = meta_output[2]
    protectdict = meta_output[3]
    retreatdict = meta_output[4]

    segmap = xsc[6]

    vars = [k for k in keys(vardict)]

    for v in vars # Iterate variables
        d = m[parse(name), parse(v)]
        
        level = vardict[v][1]
        costtype = vardict[v][2]
        println(level,costtype)

        s_arr = [segmap[i] for i in 1:size(d,1)]    # List of segments
        
        if sumsegs
            func = x -> sum(x)
            n = 1
            s_arr = ["$(size(d,1))segs"]
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
                        pval = "$(level)$(p)"
                        lev_arr = repeat([pval], outer = n)

                        outarr= [rcp_arr lev_arr s_arr cost_arr t v_arr]
                            # Write results to csv
                        if !(isfile(joinpath(outputdir,outfile)))
                            open(joinpath(outputdir,outfile),"w") do g 
                                write(g, "$(header)\n")
                            end
                        end    
             
                        open(joinpath(outputdir,outfile),"a") do g 
                            writecsv(g, outarr)
                        end

                    end
                elseif level=="retreat"
                    for j in 1:5
                        q = retreatdict[j]
                        rval = "$(level)$(q)"
                        lev_arr = repeat([rval], outer = n)

                        outarr= [rcp_arr lev_arr s_arr cost_arr t v_arr]
                        if !(isfile(joinpath(outputdir,outfile)))
                            open(joinpath(outputdir,outfile),"w") do g 
                                write(g, "$(header)\n")
                            end
                        end    
            
                        open(joinpath(outputdir,outfile),"a") do g 
                            writecsv(g, outarr)
                        end
                    end

                else
                    lev_arr = repeat([level], outer = n)
                    outarr= [rcp_arr lev_arr s_arr cost_arr t v_arr]

                    if !(isfile(joinpath(outputdir,outfile)))
                        open(joinpath(outputdir,outfile),"w") do g 
                            write(g, "$(header)\n")
                        end
                    end    
         
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
                        pval = "$(level)$(p)"
                        lev_arr = repeat([pval], outer = n)

                    elseif level=="retreat"
                        q = retreatdict[j]
                        rval = "$(level)$(q)"
                        lev_arr = repeat([rval], outer = n)
                    else
                        lev_arr = repeat([level], outer = n)
                    end

                    outarr= [rcp_arr lev_arr s_arr cost_arr t v_arr]
                    
                    if !(isfile(joinpath(outputdir,outfile)))
                        open(joinpath(outputdir,outfile),"w") do g 
                            write(g, "$(header)\n")
                        end
                    end    
                 
                    open(joinpath(outputdir,outfile),"a") do g 
                        writecsv(g, outarr)
                    end

                end
            end

        end
    end

end

function load_meta(truncate)
    metadir = "../data/meta"

    # Read header and mappings
    f = open(joinpath(metadir,"header.txt"))  # TODO define _HEADER_ path or similar a la GCAM data system
    header = readstring(f)
    close(f)

    if truncate==true
        vn = "variablenamestruncated.csv"
        println("ok")
    else
        vn = "variablenames.csv"
    end

    f = open(joinpath(metadir,vn))
    varnames = readlines(f)
    close(f)
    vardict = Dict{Any,Any}(split(varnames[i],',')[1] => (split(varnames[i],',')[2],split(varnames[i],',')[3]) for i in 1:length(varnames))

    f = open(joinpath(metadir, "protectlevelmapping.csv"))
    protect = readlines(f)
    close(f)
    protectdict = Dict{Any,Any}(parse(Int,split(protect[i],',')[1]) => parse(Int,split(protect[i],',')[2]) for i in 1:length(protect))

    f = open(joinpath(metadir, "retreatlevelmapping.csv")) 
    retreat = readlines(f)
    close(f)
    retreatdict = Dict{Any,Any}(parse(Int,split(retreat[i],',')[1]) => parse(Int,split(retreat[i],',')[2]) for i in 1:length(retreat))

    return(header,vardict, protectdict, retreatdict)
end


# Function to run tests for model
# WIP - hardcoded strings 
# Output - 
# datadir - directory where data is located; resultsdir - output directory
# gamsdata,jldata - location of comparison data from GAMS/Julia (relative to datadir)
# rcps - vector of string rcp values to test
# variables - list of variables to compare and output (strings)
function run_tests(datadir, gamsfile, jlfile, resultsdir, lslfile, subset, rcp, model=false)
    # Import model data
    modelparams = import_model_data(datadir, lslfile,"xsc.csv", subset)
    params = modelparams[1]
    xsc = modelparams[2]

    # Import GAMS Data
    gamsdata = import_comparison_data(datadir, gamsfile)
    
    # Run model for desired rcp and output results
    # TODO
    # if julia results file exists already will append to it; need to distinguish
    modellist = []


        # Run model if specified
        if model!= false
            m = run_model(params, xsc, model)
            push!(modellist,m)
        end

        # Compare and ouptut results
        jldata = import_comparison_data(datadir, jlfile) # TODO - distinct outputs for each model run

        cases = ["retreat1","retreat10","retreat100","retreat1000","retreat10000","noAdaptation","protect10",
                    "protect100","protect1000","protect10000","optimalfixed"]                     # Todo variable-length
        variables = ["protection","inundation","relocation","storms","total"]    # Todo hardcoded 
        

        for j in 1:length(cases)
            plotlist = []
            if cases[j]=="optimalfixed"
                metadata=[rcp[i],"optimalfixed","Philippines10615","total"]
                A = jldata[ (jldata[:level] .== cases[j]) .& (jldata[:costtype] .== "total"), :value] / 10 * .001
                B = gamsdata[ (gamsdata[:level] .== cases[j]) .& (gamsdata[:costtype] .== "total"), :value][2:end]
                compare_outputs(A,B,metadata,joinpath(resultsdir,"comparison.csv"))

                # Make plots
                title = string("Optimal fixed total")
                p = line_plot(B, A, title)
                savefig(p,"optimalfixed.pdf")

            else
                for k in 1:length(variables)

                    # Skip incompatible combinations
                    if (contains(cases[j], "retreat")|| cases[j]=="noAdaptation") && variables[k]=="protection"
                        continue
                    elseif contains(cases[j], "protect") && (variables[k]=="inundation" || variables[k]=="relocation")
                        continue
                    else
                        metadata=[rcp[i],cases[j],"Philippines10615",variables[k]]
                        A = jldata[ (jldata[:level] .== cases[j]) .& (jldata[:costtype] .== variables[k]), :value] / 10 * .001
                        B = gamsdata[ (gamsdata[:level] .== cases[j]) .& (gamsdata[:costtype] .== variables[k]), :value][2:end]
                        compare_outputs(A,B,metadata,joinpath(resultsdir,"comparison.csv"))
    
                        # Make plots
                        title = string(cases[j], " ", variables[k])
                        p = line_plot(B, A, title)
                        push!(plotlist, p)
                    end
                
                 end
                 make_plots(plotlist,cases[j])
            end



    end
  #  summary_report(resultsdir, "comparison.csv", resultsdir)

    if model != false
          return modellist 
    else 
        return nothing
    end

end

function run_and_write_model(datadir, resultsdir, lslfile, subset, rcp, truncate, sum)
    modelparams = import_model_data(datadir, lslfile,"xsc.csv", subset)
    params = modelparams[1]
    xsc = modelparams[2]

    m = run_model(params, xsc, "ciam")
    write_results(m, rcp, resultsdir, xsc,truncate,sum, "ciam", "results$(rcp).csv")
    return m

end

function model_driver(datadir, lslfile, subset)
    modelparams = import_model_data(datadir, lslfile,"xsc.csv", subset)
    params = modelparams[1]
    xsc = modelparams[2]

    m = run_model(params, xsc, "ciam")
    return m

end