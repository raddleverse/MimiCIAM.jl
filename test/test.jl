# Tests For CIAM Model Comparison
# TODO update comments

#------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------
using DataFrames
#using Plots
#using StatPlots
using Test
include("../src/ciamhelper.jl")



@testset "MimiCIAM" begin
    m = MimiCIAM.get_model()
    run(m)

    isneg(x)=length(x[isless.(x,0)])
    
    vargroup2D = [:WetlandNoAdapt,:FloodNoAdapt,:StormCapitalNoAdapt,:StormPopNoAdapt,:RelocateNoAdapt,
                :NoAdaptCost,:OptimalCost,:OptimalStormCapital,:OptimalStormPop,:OptimalConstruct,
                :OptimalWetland,:OptimalFlood,:OptimalRelocate,:WetlandRetreat,:WetlandProtect]
    vargroup3D = [:Construct,:StormCapitalProtect,:StormPopProtect,:StormCapitalRetreat,
                :StormPopRetreat,:FloodRetreat,:RelocateRetreat,:RetreatCost,:ProtectCost]
    vargroup4=[:OptimalH,:OptimalR,:DryLandLossOptimal]

    for v in vargroup2D
        var = m[:slrcost,v]
        @test isneg(var)==0
                    
    end

    for v in vargroup4
        var=m[:slrcost,v]
        tsteps=m[:slrcost,:ntsteps]
        for i in 2:tsteps
            @test sum(var[i,:].< var[i-1,:])==0

        end

    end
            
    for v in vargroup3D
        var=m[:slrcost,v]
        levels= size(var)[3]
        for k in 1:levels
            var2 = var[:,:,k]
            @test isneg(var2)==0
        end
    end

end






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

# Function to import comparison data and store in DataFrame
# WIP - to do -- unique header or detect header automatically
function load_data(file)
    data = readtable(file)
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


# Function to run tests for model and make plots
# WIP - hardcoded strings 
# datadir - directory where data is located; 
# resultsdir - output directory
# gamsdata,jldata - location of comparison data from GAMS/Julia (relative to datadir)
# rcp - string rcp value to test
function run_tests(gamsfile, lslfile, subset, rcp, tag, model=true, sum = false)
    # Import model data
    modelparams = import_model_data(lslfile,subset)
    params = modelparams[1]
    xsc = modelparams[2]

    # Import GAMS Data
    gamsdata = load_data(gamsfile)

    # Load Metadata
    metavars = load_meta()
    
    # Run model if specified
    if model

        m = run_model(params, xsc)

     #   if (length(subset)>5 || subset==false)
     #       sum = true
     #   else
     #       sum = false
      #  end

        write_ciam(m, xsc, rcp=rcp, tag=tag, sumsegs=sum)
    end

    # Compare and ouptut results
    jldata = load_data("../output/results-jl/results_$(rcp)_$(tag).csv") 
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
                     title = string(variables[k])
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
        make_plots(plotlist,"$(rcp)_$(levels[j])_$(tag)", resultsdir)
    end

    if model
          return m 
    else 
        return nothing
    end

end



