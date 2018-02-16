## Catherine Ledna
## Feb 9, 2018
## Model Comparison: CIAM-jl vs CIAM-GAMS
## Single segment -- Philippines

using Mimi
using GR
using Plots
include("../src/ciam_debug.jl")
include("../src/helper.jl")

# Step 1: Import and run RCP8.5 single-segment scenario in ciam_debug
data_dir = joinpath("test_phil/input-data")
maindata =  ["philinput.csv", "globalparams.csv","pop.csv","ypcc.csv","ypc_usa.csv"]
lsldata = "philinputlsl_reshape.csv"
  
params = Dict{Any, Any}(lowercase(splitext(m)[1]) => readdlm(joinpath(data_dir,m), ',' ) for m in maindata)
lslparams = Dict{Any, Any}(lowercase(splitext(lsldata)[1]) => readdlm(joinpath(data_dir,lsldata), ',' ))

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
lslrall = Dict{Any,Any}()
parse_long!(lslparams["philinputlsl_reshape"], lslrall,1,1, true)



mainparams = Dict{Any,Any}()
parse_long!(params["philinput"],mainparams,1)
parse_long!(params["globalparams"],mainparams)
mainparams["ypc_usa"] = params["ypc_usa"][2:21,2]
mainparams["pop"] = params["pop"][2:21,2:end]
mainparams["ypcc"] = params["ypcc"][2:21,2:end]
mainparams["lslr"] = lslrall["rcp85_p50"]
a=["refpopdens","popdens","gtapland","length","cci","psig0","psig0coef","psigA","psigB","rsig0","rsigA",
    "rsigB","wetland","s10","s100","s1000","smax","area1","area2","area3","area4","area5","area6","area7","area8",
    "area9","area10","area11","area12","area13","area14","area15"] 
for k in a
    mainparams[k] = [convert(Float64,mainparams[k])]
end

b = ["movefactor","depr","kgdp"]
for k in b
    mainparams[k] = convert(Float64,mainparams[k])
end

#mainparams["wetland"] = [0.0]
 
xsc = "xsc.csv"
xsc_params = Dict{Any, Any}(lowercase(splitext(xsc)[1]) => readdlm(joinpath(data_dir,xsc), '\r' ))
xscout = prepxsc(xsc_params)
xsc_ind = xscout[1]
allrgns = sort(xscout[2])
allsegs = sort(xscout[3])


# Run model 
function run_model(params, xsc)
    atpers = [1, 5, 10, 15, 19]
    m = Model()
    setindex(m, :time, 20)
    setindex(m, :adaptPers, length(atpers))
    setindex(m, :regions, allrgns)
    setindex(m, :segments, allsegs)

    addcomponent(m, ciam)

    setparameter(m, :ciam, :ntsteps, 20)       
    setparameter(m, :ciam, :tstep, 10.)
    setparameter(m, :ciam, :at, atpers) 

    adapt = [1,10,100,1000,10000]
    setparameter(m, :ciam, :adaptOptions, adapt)
    setparameter(m, :ciam, :fixed, true)
    setparameter(m, :ciam, :landinput, false)

    setparameter(m, :ciam, :xsc, xsc)
    setleftoverparameters(m, params)

    run(m)

    return m
end

n = run_model(mainparams, xsc_ind)

function write_results(model, dirname, outdirname, flood=false)
    dir = joinpath(outdirname, dirname)

    if isdir(dir)
    else
        mkdir(dir)
    end

    writecsv(joinpath(dir, "RegionalCost.csv"), model[:ciam, :RegionalCost])
    writecsv(joinpath(dir, "SegCost.csv"), model[:ciam, :AdaptationCost])
    writecsv(joinpath(dir, "NoAdaptCost.csv"), model[:ciam, :NoAdaptCost])
    writecsv(joinpath(dir, "WetlandNoAdaptCost.csv"), model[:ciam, :WetlandNoAdapt])
    writecsv(joinpath(dir, "FloodNoAdaptCost.csv"), model[:ciam, :FloodNoAdapt])
    writecsv(joinpath(dir, "StormNoAdaptCost.csv"), model[:ciam, :StormNoAdapt])
    writecsv(joinpath(dir, "RelocateNoAdaptCost.csv"), model[:ciam, :RelocateNoAdapt])

    if flood==true
        writecsv(joinpath(dir, "landvalue.csv"), model[:ciam, :landvalue])
        writecsv(joinpath(dir, "coastarea.csv"), model[:ciam, :coastArea])
        writecsv(joinpath(dir, "capital.csv"), model[:ciam, :capital])
        writecsv(joinpath(dir, "ypc.csv"), model[:ciam, :ypc_seg])
        writecsv(joinpath(dir, "popdens.csv"), model[:ciam, :popdens])
    
    
    end

    
end

write_results(n, "rcp85_p50")

# Run for rcp0_p5, rcp0_p50, rcp0_p95
# Updated round 2: Fixed area function
rcps = ["rcp0_p5", "rcp0_p50", "rcp0_p95", "rcp85_p50"]
for r in 1:4
    paramdir = copy(mainparams)
    paramdir["lslr"] = lslrall[rcps[r]]

    out = run_model(paramdir, xsc_ind)
    write_results(out, rcps[r], "test_phil/results-jl/round2",true)
end



# Step 2: Compare outputs to GAMS version

