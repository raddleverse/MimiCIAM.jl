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
    adaptperiods = [1, 5, 10, 15, 19]
    m = Model()
    setindex(m, :time, 20)
    setindex(m, :adaptPers, length(adaptperiods))
    setindex(m, :regions, xsc[2])
    setindex(m, :segments, xsc[3])

    addcomponent(m, eval(parse(model)))

    setparameter(m, parse(model), :ntsteps, 20)       
    #setparameter(m, :ciam, :tstep, 10.)
    setparameter(m, parse(model), :at, adaptperiods) 

    adaptoptions = [1,10,100,1000,10000]
    setparameter(m, parse(model), :adaptOptions, adaptoptions)
    setparameter(m, parse(model), :fixed, true)
    setparameter(m, parse(model), :landinput, true)

    setparameter(m, parse(model), :xsc, xsc[1])
    setleftoverparameters(m, params)

    run(m)

    return m
end

# Helper function to trim arrays based on specifications
# Used for importing data
# data - input data; outdict - output dictionary that is modified by fcn
# csnip = # of columns to remove from beginning; rsnip = # of rows to remove from beginning
# arr = bool for whether to output as 2-d array or 1-d vector 
# TODO: if we standardize data preferences in main CIAM, arr is prob less necessary?
function parse_long!(data::Array{Any,2}, outdict::Dict{Any,Any}, csnip = 0, rsnip=0,arr=false)
    data = data[(rsnip+1):end, (csnip+1):end ]
    
    for i in 1:size(data,1)
        if arr
            outdict[data[i,1]] = data[i:i,2:end]
        else
            outdict[data[i,1]] = data[i,2:end][1] 
        end
        
    end

end


# Wrapper for importing model data. Not generalizable; hard-coded names; WIP
# datadir - data directory
# paramfiles - filenames for main parameters (string)
# lslfile - filename for lsl (string)
# xscfile - filename for segment-country mapping (string)
function import_model_data(datadir, paramfiles, lslfile, xscfile)
    # Process main and lsl params
    params = Dict{Any, Any}(lowercase(splitext(m)[1]) => readdlm(joinpath(datadir,m), ',' ) for m in paramfiles)
    lslparams = Dict{Any, Any}(lowercase(splitext(lslfile)[1]) => readdlm(joinpath(data_dir,lslfile), ',' ))
    

    # Main params 
    mainparams = Dict{Any,Any}()
    parse_long!(params["philinput"],mainparams,1)
    parse_long!(params["globalparams"],mainparams)
    mainparams["ypc_usa"] = params["ypc_usa"][2:21,2]
    mainparams["pop"] = params["pop"][2:21,2:end]
    mainparams["ypcc"] = params["ypcc"][2:21,2:end]
    
    # for p in keys(mainparams)
    #     if isinteger(mainparams[p])
    #         mainparams[p] = [convert(Float64,mainparams[p])]
    #     end
    # end
    a=["refpopdens","popdens","gtapland","length","cci","psig0","psig0coef","psigA","psigB","rsig0","rsigA",
    "rsigB","wetland","s10","s100","s1000","smax","area1","area2","area3","area4","area5","area6","area7","area8",
    "area9","area10","area11","area12","area13","area14","area15"] 
    for k in a
        mainparams[k] = [convert(Float64,mainparams[k])]
    end
    b = ["movefactor","depr","kgdp","tstep"]
    for k in b
        mainparams[k] = convert(Float64,mainparams[k])
    end


    # LSLR Dictionary
    lslrall = Dict{Any,Any}()
    parse_long!(lslparams["philinputlsl_reshape"], lslrall,1,1, true)

    # Process XSC
    xsc_params = Dict{Any, Any}(lowercase(splitext(xscfile)[1]) => readdlm(joinpath(datadir,xscfile), '\r' ))
    xsc = prepxsc(xsc_params)

    return(mainparams, xsc, lslrall)

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

# Function to run tests for model
# WIP - hardcoded strings 
# Output - 
# datadir - directory where data is located; resultsdir - output directory
# gamsdata,jldata - location of comparison data from GAMS/Julia (relative to datadir)
# rcps - vector of string rcp values to test
# variables - list of variables to compare and output (strings)
function run_tests(datadir, paramfiles, gamsfile, jlfile, resultsdir, rcps, model=false)
    # Import model data
    modelparams = import_model_data(datadir, paramfiles, "philinputlsl_reshape.csv","xsc.csv")
    params = modelparams[1]
    xsc = modelparams[2]
    lslrall = modelparams[3]

    # Import GAMS Data
    gamsdata = import_comparison_data(datadir, gamsfile)
    
    # Run model for desired RCPs and output results
    # TODO 
    # if julia results file exists already will append to it; need to distinguish
    modellist = []
    for i in 1:length(rcps)

        # Run model if specified
        if model!= false
            params["lslr"] = lslrall[rcps[i]]
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
                metadata=[rcps[i],"optimalfixed","Philippines10615","total"]
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
                        metadata=[rcps[i],cases[j],"Philippines10615",variables[k]]
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


    end
  #  summary_report(resultsdir, "comparison.csv", resultsdir)

    if model != false
        println("returning model")
        return modellist 
    else 
        return nothing
    end

end

