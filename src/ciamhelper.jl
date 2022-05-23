using DataFrames
using CSV
using Query

## ---------------------------------
## 1. Read In Helper Functions
## ---------------------------------

"""
    load_ciam_params()

Load CIAM parameters from CSV to dictionary.
"""
function load_ciam_params()

    data_dir = joinpath(@__DIR__, "..","data","input")
    files = readdir(data_dir)
    filter!(i->(i!="desktop.ini" && i!=".DS_Store" && i!="xsc.csv"), files)

    params = Dict{Any, Any}(lowercase(splitext(m)[1]) => CSV.read(joinpath(data_dir,m),DataFrame) |> DataFrame for m in files)

    return params

end

"""
    preplsl!(lslfile,subset, params,segnames)

Read in LSLR from file and filter to desired set of segments, note that this modifies
the input parameter dictionary `params`.  The arguments are as fullows:

- lslfile - name of lslr file to use; location relative to data/input-data directory
- subset - list of segments you want
- params - parameter dictionary you want to add lslr to
- segnames - names of segments
"""
function preplsl!(lslfile, subset, params, segnames)

    data_dir = joinpath(@__DIR__,"..","data","lslr")
    lsl_params = CSV.read(joinpath(data_dir,lslfile), DataFrame) |> DataFrame

    # Filter according to subset segments
    if subset != false
        col_names = [i for i in names(lsl_params) if string(i) in subset]
        lsl_params = lsl_params[!,col_names]
    end

    # Chomp off unrelated rows and sort alphabetically (do this regardless of whether there's a subset)
    col_names = [i  for i in names(lsl_params) if string(i) in segnames]
    col_names = sort(col_names)
    lsl_params = lsl_params[!,col_names]

    params["lslr"] = convert(Array{Float64,2}, Matrix(lsl_params))

    return params
end

"""
    prepssp!(ssp, ssp_simplified, params, rgnnames, segnames, popinput)

Read in SSP from file and filter to desired set of segments, note that this modifies
the input parameter dictionary `params`.  The arguments are as fullows:

- ssp - name of lslr file to use; location relative to data/input-data directory
- ssp_simplified - list of segments you want
- params - parameter dictionary you want to add lslr to
- rgnnames - names of regions (ciam_country)
- segnames - names of segments
- popinput - population density data set (0=original CIAM, 1=Jones and O'Neill 2016 (not supported), 2=Merkens et al 2016 (not supported))
"""
function prepssp!(ssp, ssp_simplified, params, rgnnames, segnames, popinput)

    data_dir = joinpath(@__DIR__,"..","data","ssp")

    # read and set population densities for Jones and Merkens data sets whether or not
    # we are using  them, so they have some defaults.
    # popinput=0 is the only supported option at this time
    if popinput != 0

        error("The `popinput` argument values of 1 and 2 are not supported at this time.  In the future they will indicate use of Jones and O'Neill 2016 or Merkens et al 2016 population data, respectively.")

        popdens_seg_jones = CSV.read(joinpath(data_dir,string("popdens_seg_jones_ssp",ssp_simplified,".csv")), DataFrame)
        popdens_seg_merkens=CSV.read(joinpath(data_dir,string("popdens_seg_merkens_ssp",ssp_simplified,".csv")), DataFrame)

        seg_col_names = [i for i in names(popdens_seg_jones) if string(i) in segnames]
        sort!(seg_col_names)
        popdens_seg_jones = popdens_seg_jones[!, seg_col_names]

        seg_col_names = [i for i in names(popdens_seg_merkens) if string(i) in segnames]
        sort!(seg_col_names)
        popdens_seg_merkens = popdens_seg_merkens[!, seg_col_names]

        params["popdens_seg_jones"]     = Array{Float64,2}(popdens_seg_jones)
        params["popdens_seg_merkens"]   = Array{Float64,2}(popdens_seg_merkens)
    end

    if ssp == false # Do nothing, base ssp data already loaded
        return params

    else
        pop = CSV.read(joinpath(data_dir,string("pop_",ssp,".csv")), DataFrame)
        ypc = CSV.read(joinpath(data_dir, string("ypcc_",ssp,".csv")), DataFrame)

        col_names = [i  for i in names(pop) if string(i) in rgnnames]
        sort!(col_names)
        pop = pop[!, col_names]
        ypc = ypc[!, col_names]

        params["pop"]   = Array{Float64,2}(pop)
        params["ypcc"]  = Array{Float64,2}(ypc)
    end
    return params
end

"""
    parse_ciam_params!(params, rgn_order, seg_order, surgeoption)

Process CIAM data from csv to usable format and store outputs in params, note that
this modifies the input parameter dictionary `params`, and the funciton is
specific to CIAM data so there are some hard-coded names and assumptions. The
arguments are as fullows:

- rgn_order - alphabetized lists of regions used (ciam_country)
- seg_order - alphabetized lists of segments used
- surgeoption - which surge exposure data set to be used
"""
function parse_ciam_params!(params, rgn_order, seg_order, surgeoption)

    # we need to grab the original keys so it doesn't try to recurse when we
    # make new entries into the dictionary
    original_keys = [k for k in keys(params)]

    for k in original_keys
        p = params[k] # Data frame

        # data key case
        if k == "data"
            colnames = filter(f -> string(f) != "NA", names(p)) # Preserve column names

            # Filter segments to subset
            segs = p[!,1]
            seg_inds = filter_index(segs, seg_order)
            p = p[seg_inds,:]

            # Sort alphabetically
            seg_alpha = sortperm(p[!,1])
            p = p[seg_alpha,:]

            # Process all variables
            if length(seg_inds) >= 1
                for k in 2:(length(colnames) + 1)
                    varname = string(colnames[k-1])
                    newvars = p[1:end,k]
                    newvars = [convert(Float64,v) for v in newvars]
                    params[varname] = newvars
                end
                delete!(params, "data")
            else
                error("Segments in dictionary do not match supplied segments")
            end

        # globalparams key case
        elseif k == "globalparams"

            for k in 1:nrow(p)

                varname = p[k, 1]
                newval = p[k, 2]

                if (varname == "ntsteps" || varname == "adaptPers")
                    newval = parse(Int64,newval)
                elseif lowercase(newval) == "true"
                    newval = true
                elseif lowercase(newval) == "false"
                    newval = false
                else
                    newval = parse(Float64,newval)
                end
                params[varname] = newval
            end
            delete!(params, "globalparams")

        # surge exposure key case
        elseif occursin("surgeexposure",k) # generalize for multiple possible surge data sets

            p = @from i in p begin
                @where i.segments in seg_order
                @select i
                @collect DataFrame
            end

            # Sort alphabetically
            sort!(p, :segments)
            params[k] = convert(Array{Float64,2}, Matrix(p[:,2:6]))

        # refa key case
        elseif k == "refa_h" || k == "refa_r"

            p = @from i in p begin
                @where i.segments in seg_order
                @select i
                @collect DataFrame
            end

            sort!(p, :segments)
            delete!(params, k)

            if k=="refa_h"
                params["refA_H"] = convert(Array{Float64}, p[!,:value])
            else
                params["refA_R"] = convert(Array{Float64}, p[!,:value])
            end

        # Country data matrices parameter case
        elseif size(p, 2) == 2

            # Filter regions
            r_inds = filter_index(p[!,1], rgn_order)
            p = p[r_inds,:]

            # Alphabetize
            p = p[sortperm(p[!,1]),:]
            if p[!,1] != rgn_order
                error("Regions in dictionary do not match supplied regions, ", k)
            else
                newvals = p[!,2]
                params[k] = Array{Float64,1}(newvals) # Coerce to Array{Float64,1}
            end

        # Time-country data matrices parameter case
        elseif size(p, 2) > 3

            # Alphabetize
            col_names = [i for i in names(p) if string(i) in rgn_order]
            p = p[!, sort(col_names)]
            params[k] = Array{Float64,2}(p)

        # Single dimension parameter case
        elseif size(p, 2)==1
            params[k] = Array{Float64,1}(p[!,1])
        end
    end

    # set p["surgeexposure"] to the one designated by the surgeoption argument
    # by default the original data file is read in as p["surgeexposure"], so in
    # that case no action needed. Otherwise...
    if surgeoption == 1
        params["surgeexposure"] = params["surgeexposure_dc-gtsr"]
    elseif surgeoption == 2
        params["surgeexposure"] = params["surgeexposure_gtsr"]
    elseif surgeoption != 0
        error("The `surgeoption` argument must be 0, 1 or 2.")
    end

end

"""
    filter_index(v1, v2)

Filter a vector (v1) by a second vector (v2) and return indices of contained elements
"""
function filter_index(v1, v2)
    out = []
    for i in 1:length(v1)
        if v1[i] in v2
            push!(out, i)
        end
    end
    return(out)
end

"""
Process the segment-country mapping file (xsc) in CIAM by (1) Reads from CSV
and outputs list of dictionaries and arrays (2) Filters xsc file to desired
segments/regions
"""
function prepxsc(subset)

    data_dir = joinpath(@__DIR__,"..","data","input")
    xscfile = "xsc.csv"
    xsc_name = replace(xscfile, r".csv" => s"") # Strip ".csv" from file name

    # Read in csv and convert to dictionary format
    xsc_params = Dict{Any, Any}(lowercase(splitext(xscfile)[1]) => CSV.read(joinpath(data_dir, xscfile), DataFrame))
    xsc_char = Dict{Any,Any}( xsc_params[xsc_name][i,1] => (xsc_params[xsc_name][i,2],xsc_params[xsc_name][i,3], xsc_params[xsc_name][i,4]) for i in 1:size(xsc_params[xsc_name],1))

    # If only a subset of segments is used, filter down to relevant segments
    if subset!=false
        filter!(p-> (first(p) in subset), xsc_char)
    end

    # Create region and segment indices
    rgns = sort(unique([i[1] for i in collect(values(xsc_char))]))
    segs = string.(sort(unique(collect(keys(xsc_char)))))

    xsc_ind = Dict{Any,Any}()      # numeric seg -> (numeric rgn, greenland bool)
    xsc_segmap = Dict{Any,Any}()   # Numeric seg/rgn -> char seg/rgn
    xsc_rgnmap = Dict{Any,Any}()

    for i in 1:length(segs)
        r = xsc_char[segs[i]][1]   # Region character
        grn = xsc_char[segs[i]][2] # 0 = non-Greenland, 1 = greenland bool
        isl = xsc_char[segs[i]][3] # 0 = non-island, 1 = island bool
        r_ind = findind(r, rgns)   # Region index

        new_val = (r_ind, grn, isl)     # New tuple w/ region index instead of character

        # Build XSC Seg/rgn Maps
        r2 = rgns[r_ind]           # New region char
        s = segs[i]
        xsc_segmap[i] = s
        if !(r2 in values(xsc_rgnmap))
            xsc_rgnmap[r_ind] = r2
        end

        xsc_ind[i] = new_val
    end

    return (xsc_ind, rgns, segs, xsc_rgnmap)

end

"""
    findind(val, vec)
Look up index corresponding to name with arguments `vec`, a vector of region or
segment names (strings) and `val`, a string corresponding to value in 'vec'
"""
function findind(val, vec)
    name_ind = findfirst(isequal(val), vec)
    return name_ind
end

function load_subset(subset=false)
    dir = joinpath(@__DIR__,"..","data","subsets")
    if subset != false
        subs = readlines(joinpath(dir,subset))
        return subs
    else
        return false
    end
end

function init(; f::Union{String, Nothing} = nothing)
    if isnothing(f)
        f = joinpath(@__DIR__,"..","data","batch","init.csv")
    end
    varnames=CSV.read(f, DataFrame)
    vardict = Dict{Any,Any}( String(i) => varnames[!,i] for i in names(varnames))
    return(vardict)
end

"""
    import_model_data(lslfile, sub, ssp, ssp_simplified, popinput)

Wrapper for importing model data with the arguments:
- lslfile - filename for lsl (string)
- sub - filename with names of segments to use (string) or false (bool) to run all segments
- ssp - SSP scenario + modeling group (specific long name)
- ssp_simplified  - SSP scenario (1-5)
- popinput - population density data set to be used (0=original CIAM, 1=Jones and O'Neill 2016 (not supported), 2=Merkens et al 2016 (not supported))
- surgeoption - surge exposure data set to use (0=original CIAM/DINAS-COAST, 1=D-C corrected by GTSR/D-C bias, 2=GTSR nearest data point(s))
"""
function import_model_data(lslfile,sub,ssp,ssp_simplified,popinput,surgeoption)

    subset = (sub == "false") ? false : load_subset(sub)

    # Process main and lsl params
    params = load_ciam_params()

    # Process XSC (segment-country mapping dictionary)
    xsc = prepxsc(subset)

    # Process params using xsc and format lsl file
    parse_ciam_params!(params, xsc[2], xsc[3], surgeoption)
    preplsl!(lslfile, subset, params,xsc[3])
    prepssp!(ssp,ssp_simplified,params,xsc[2],xsc[3],popinput)

    return(params, xsc)

end

function load_meta()
    metadir = joinpath(@__DIR__,"..","data","meta")

    # Read header and mappings
    header = read(open(joinpath(metadir,"header.txt")),String)

    varnames = readlines(open(joinpath(metadir,"variablenames.csv")))
    vardict = Dict{Any,Any}(split(varnames[i],',')[1] => (split(varnames[i],',')[2],split(varnames[i],',')[3]) for i in 1:length(varnames))

    protect = readlines(open(joinpath(metadir, "protectlevels.csv")))
    protectdict = Dict{Any,Any}(parse(Int,split(protect[i],',')[1]) => split(protect[i],',')[2] for i in 1:length(protect))

    retreat = readlines(open(joinpath(metadir, "retreatlevels.csv")))
    retreatdict = Dict{Any,Any}(parse(Int,split(retreat[i],',')[1]) => split(retreat[i],',')[2] for i in 1:length(retreat))

    return(header,vardict, protectdict, retreatdict)
end


## ---------------------------------
## 2. Write Out Helper Functions
## ---------------------------------

"""
    write_ciam(m; runname::String = "base", sumsegs::String = "seg", varnames::Bool = false, tag::Bool = false)

Write out model results to CSV file using arguments:

- model - output from get_model()
- runname
- sumsegs - whether to sum across all segments, to region level, or no sums
- varnames - if not false, write the passed variable names; if false get defaults from file
- tag
To do: possibly modify to work with DataVoyager()
"""
function write_ciam(model; outputdir::String = joinpath(@__DIR__,"..","output"),runname::String = "base", sumsegs::String = "seg", varnames::Bool = false, tag::Bool = false)

    meta_output = load_meta()
    rcp         = model[:slrcost,:rcp]
    pctl        = model[:slrcost,:percentile]
    ssp         = model[:slrcost,:ssp]
    fixed       = (model[:slrcost, :fixed]) ? "fixed" : "flex"
    rcp_str     = "$(rcp)p$(pctl)ssp$(ssp)$(fixed)"

    xsc         = load_xsc()
    segmap      = load_segmap()
    segRgnDict  = Dict{Any,Any}(xsc[!,:seg][i] => xsc[!,:rgn][i] for i in 1:size(xsc,1))

    if varnames == false
        varnames = [k for k in keys(meta_output[2])] # to do change
    end

    vargroup1 = [] # 2D vars
    vargroup2 = [] # vars greater than 2D

    for v in varnames
        if length(size(model[:slrcost,Symbol(v)])) > 2
            push!(vargroup2, Symbol(v))
        else
            push!(vargroup1, Symbol(v))
        end
    end

    # Assign 2D variables to dataframe
    # 2 cases: 1. adapt pers is first; 2. adapt pers is second
    common_order = [:time, :ciam_country, :segments, :level]

    for i in 1:length(vargroup1)
        temp = getdataframe(model, :slrcost, vargroup1[i])
        missing_names = [j for j in common_order if !(String(j) in names(temp))]

        if length(missing_names) >= 1
            for name in missing_names
                temp[!,name] = missings(size(temp)[1])
            end
        end

        if :ciam_country in missing_names && !(:segments in missing_names)
            temp = temp |> @map(merge(_,{ciam_country=segRgnDict[_.segments]})) |> DataFrame
        end

        temp[!,:variable] = fill(String(vargroup1[i]), nrow(temp))
        rename!(temp,vargroup1[i] => :value)
        temp = temp[!, [:time, :ciam_country, :segments, :level, :variable, :value]]

        if i == 1
            global df = temp
        else
            df = [df; temp]
        end
    end

    # Assign 3D variables to second data frame and join
    ntime = model[:slrcost, :ntsteps]
    segID = model[:slrcost, :segID]
    colnames = Symbol.(segID_to_seg(Int64.(segID),segmap))

    for j in 1:length(vargroup2)

        ndim1 = size(model[:slrcost,vargroup2[j]])[3]

        for k in 1:ndim1

            temp = DataFrame(model[:slrcost,vargroup2[j]][:,:,k], :auto)

            rename!(temp, colnames )
            temp[!,:time] = 1:ntime

            if String(vargroup2[j])=="Construct" || occursin("Protect",String(vargroup2[j]))
                dim1 = k+1
                adapt=model[:slrcost,:adaptoptions][dim1]
            else
                dim1=k
                adapt=model[:slrcost,:adaptoptions][dim1]
            end

            temp[!,:level] = fill(adapt, ntime)
            temp = stack(temp,colnames)

            rename!(temp,:variable => :segments)

            temp[!,:segments] = [String(i) for i in temp[!,:segments]]
            temp[!,:variable]= fill(String(vargroup2[j]),nrow(temp))

            temp = temp |> @map(merge(_,{ciam_country=segRgnDict[_.segments]})) |> DataFrame
            temp = temp[!,[:time,:ciam_country,:segments,:level,:variable,:value]]

            if j == 1 && k == 1
                global df2 = temp
            else
                df2 = [df2; temp]
            end
        end
    end

    # Sum to either region-level, global-level, or leave as seg-level
    outdf = [df; df2]
    outfile = tag ? "$(runname)_$(sumsegs)_$(rcp_str)_$(tag).csv" : "$(runname)_$(sumsegs)_$(rcp_str).csv"

    if sumsegs=="rgn"
        rgndf = outdf |> @groupby({_.time,_.ciam_country, _.level, _.variable}) |> @map(merge(key(_),{value = sum(_.value)})) |> DataFrame
        CSV.write(joinpath(outputdir,outfile),rgndf)
    elseif sumsegs=="seg"
        CSV.write(joinpath(outputdir,outfile),outdf)
    elseif sumsegs=="global"
        globdf = outdf |> @groupby({_.time, _.level, _.variable}) |>
            @map(merge(key(_),{value = sum(_.value)})) |> DataFrame
        CSV.write(joinpath(outputdir,outfile),globdf)
    end
end

"""
    write_optimal_costs(model; outputdir::String = joinpath(@__DIR__,"..","output"), runname="base")

Streamline writing results for optimal adaptation costs.
"""
function write_optimal_costs(model; outputdir::String = joinpath(@__DIR__,"..","output"), runname="base")

    # Output: Data Frame with segment,region,time,level,option, suboption
    #   E.g. 'OptimalProtect', 'Construct'
    # Should output 2 CSVs: 1 with just the 3 main categories, 2nd with
    #   detailed subcategories

    rcp     = model[:slrcost,:rcp]
    pctl    = model[:slrcost,:percentile]
    ssp     = model[:slrcost,:ssp]
    fixed   = (model[:slrcost, :fixed]) ? "fixed" : "flex"
    rcp_str = "$(rcp)p$(pctl)ssp$(ssp)$(fixed)"
    xsc     = load_xsc()
    segRgnDict = Dict{Any,Any}( xsc[!,:seg][i] => xsc[!,:rgn][i] for i in 1:size(xsc,1))

    # 1. Create aggregate adaptation decision DF
    temp1 = getdataframe(model, :slrcost => :OptimalCost)
    temp1 = temp1 |> @map(merge(_,{ciam_country=segRgnDict[_.segments]})) |> DataFrame

    temp2 = getdataframe(model, :slrcost => :OptimalLevel)
    temp3 = getdataframe(model, :slrcost => :OptimalOption)

    # Join dataframes and reorganize
    # Inner join to restrict to only rows in both dataframes
    # (which should be all of time, because all segments should have values at all time steps)
    out = innerjoin(temp1,temp2, on = [:time,:segments])
    out = innerjoin(out,temp3, on = [:time,:segments])

    # Replace OptimalOption numeric value with string
    lookup = Dict{Any,Any}(-2.0=> "RetreatCost", -1.0=> "ProtectCost",-3.0=>"NoAdaptCost")
    out = out |> @map(merge(_,{variable=lookup[_.OptimalOption]})) |> DataFrame
    rename!(out, Dict(:OptimalLevel => :level))
    out = out[!,[:time, :ciam_country, :segments, :variable, :level, :OptimalCost]]

    # Write to file
    outfile = "$(runname)_seg_$(rcp_str)_optimal.csv"
    CSV.write(joinpath(outputdir,outfile),out)

    # Write Sub-Costs
    vars = [:OptimalStormCapital, :OptimalStormPop, :OptimalConstruct,
            :OptimalFlood, :OptimalRelocate, :OptimalWetland]

    for i in 1:length(vars)

        temp = getdataframe(model, :slrcost => vars[i])
        temp = temp |> @map(merge(_,{ciam_country=segRgnDict[_.segments]})) |> DataFrame

        temp[!,:variable]= fill(String(vars[i]),nrow(temp))

        temp2 = getdataframe(model, :slrcost => :OptimalLevel)
        temp3 = getdataframe(model, :slrcost => :OptimalOption)

        # Join dataframes and reorganize
        out = innerjoin(temp,temp2, on=[:time,:segments])
        out = innerjoin(out,temp3, on=[:time,:segments])

        # Replace OptimalOption numeric value with string
        lookup = Dict{Any,Any}(-2.0=> "RetreatCost", -1.0=> "ProtectCost",-3.0=>"NoAdaptCost")
        out = out |> @map(merge(_,{AdaptCategory=lookup[_.OptimalOption]})) |> DataFrame
        rename!(out, Dict(:OptimalLevel => :level))
        rename!(out,vars[i]=>:value)
        out = out[!,[:time,:ciam_country,:segments,:AdaptCategory,:variable,:level,:value]]

        if i==1
            global df = out
        else
            df = [df;out]
        end
    end

    # Write to CSV
    outfile2 = "$(runname)_seg_$(rcp_str)_optimal_subcost.csv"
    CSV.write(joinpath(outputdir,outfile2),df)

end

function write_optimal_protect_retreat(model; runname="base")
    outputdir = joinpath(@__DIR__,"..","output")
    rcp     = model[:slrcost,:rcp]
    pctl    = model[:slrcost,:percentile]
    ssp     = model[:slrcost,:ssp]
    fixed   = (model[:slrcost,:fixed]) ? "fixed" : "flex"
    rcp_str = "$(rcp)p$(pctl)ssp$(ssp)$(fixed)"
    xsc     = load_xsc()
    # segRgnDict = Dict{Any,Any}( xsc[!,:seg][i] => xsc[!,:rgn][i] for i in 1:size(xsc,1)) # note this isn't used

    # 1. Create aggregate adaptation decision DF
    pl = getdataframe(model, :slrcost => :OptimalProtectLevel)
    rl = getdataframe(model, :slrcost => :OptimalRetreatLevel)
    out = innerjoin(pl,rl, on=[:time,:segments])

    ## To do: read in protect, retreat variables and filter to optimal levels; add to out
    protect = model[:slrcost, :ProtectCost]
    retreat = model[:slrcost, :RetreatCost]
    for i in 1:5
        if i>1
            prot = DataFrame(model[:slrcost,:ProtectCost][:,:,i-1])
            ret = DataFrame(model[:slrcost,:RetreatCost][:,:,i])
        else
            ret = DataFrame(model[:slrcost,:RetreatCost][:,:,i])
        end
    end

    outfile = "$(runname)_seg_$(rcp_str)_ProtectRetreat.csv"
    CSV.write(joinpath(outputdir,outfile),out)

end

## ---------------------------------
## 3. Basic Functions for Segment-Region Lookup
## ---------------------------------

function load_segmap()
    segmap = CSV.read(joinpath(@__DIR__, "..", "data","meta", "segIDmap.csv"), DataFrame) |> DataFrame
    segmap[!,:segID] = [Int64(i) for i in segmap[!,:segID]]
    return(segmap)
end

function load_xsc()
    xsc = CSV.read(joinpath(@__DIR__,"..","data","input","xsc.csv"), DataFrame) |> DataFrame
    return(xsc)
end

"""
    segID_to_seg(segID, segmap)

Look up string name of segment from segID and return only first result for each
ID entry as an array. The arguments are as follows:

- segID - int or array of ints
- segmap - output of load_rgnmap or load_segmap() (DataFrame)
"""
function segID_to_seg(segID, segmap)
    seg = [String(segmap[!,:seg][segmap.segID.==i][1]) for i in segID]
    return(seg)
end

"""
    segStr_to_segID(segstr)

Get segID from string name of segment. Segstr must be an array of strings
"""
function segStr_to_segID(segstr)
    ids = [parse(Int64,replace(i, r"[^0-9]"=> "")) for i in segstr]
    return(ids)
end

"""
    getTimeSeries(model, ensnum; segIDs = false, rgns = false, sumsegs = "global")

Return time series of costs at global, regional and/or segment level scale and also
compute cost as percent of regional or global gdp.
"""
function getTimeSeries(model, ensnum; segIDs = false, rgns = false, sumsegs = "global")

    # If not using segment-level aggregation, segIDs refers to
    # individual segments to report in addition to global/regional

    if sumsegs == "seg" # Report all segments in model or those specified
        if segIDs == false
            segIDs = model[:slrcost,:segID]
        end
    end

    xsc = MimiCIAM.load_xsc()
    segRgnDict = Dict{Any,Any}( xsc[!,:seg][i] => (xsc[!,:rgn][i],xsc[!,:segID][i]) for i in 1:size(xsc,1))

    # Write Main and Sub-Costs
    vars = [:OptimalCost,:OptimalStormCapital, :OptimalStormPop, :OptimalConstruct,
        :OptimalFlood, :OptimalRelocate, :OptimalWetland]
    global df = DataFrame()

    for i in 1:length(vars)

        temp = MimiCIAM.getdataframe(model, :slrcost => vars[i])
        temp = temp |> @map(merge(_,{ciam_country=segRgnDict[_.segments][1],segID=segRgnDict[_.segments][2]})) |> DataFrame
        #temp[!,:costtype]= String(vars[i])

        temp2 = MimiCIAM.getdataframe(model, :slrcost => :OptimalLevel)
        temp3 = MimiCIAM.getdataframe(model, :slrcost => :OptimalOption)
        temp4 = MimiCIAM.getdataframe(model, :slrcost => :ypcc)
        temp5 = MimiCIAM.getdataframe(model, :slrcost => :pop)

        # Join dataframes and reorganize
        out = innerjoin(temp,temp2, on=[:time,:segments])
        out = innerjoin(out,temp3, on=[:time,:segments])
        out = innerjoin(out,temp4, on=[:time,:ciam_country])
        out = innerjoin(out,temp5, on=[:time,:ciam_country])

        # Replace OptimalOption numeric value with string
        lookup = Dict{Any,Any}(-2.0=> "Retreat", -1.0=> "Protection",-3.0=>"No Adaptation")
        out = out |> @map(merge(_,{category=lookup[_.OptimalOption]})) |> DataFrame
        rename!(out, Dict(:OptimalLevel => :level))
        rename!(out,vars[i]=>:cost)

        out[!,:ens] = fill(ensnum, size(out)[1])
        col_order=[:ens,:time,:ciam_country,:segments,:segID,:category,:level,:cost,:ypcc,:pop]
        out = out[!,col_order]

        # Aggregate to geographic level
        subset = filter(row -> row[:segID] in segIDs,out)

        if sumsegs=="rgn"
            rgndf = out |> @groupby({_.ens,_.time,_.ciam_country, _.level, _.category,_.ypcc,_.pop}) |> @map(merge(key(_),{cost = sum(_.cost)})) |> DataFrame
            rgndf[!,:segID] = fill(0., size(rgndf)[1])
            rgndf[:segments] = fill("regional", size(rgndf)[1])
            rgndf = rgndf[!,col_order]
            rgndf = [rgndf;subset]

            rgndf[!,:gdp]=rgndf[!,:ypcc].*rgndf[!,:pop]./1e3 # GDP is in $Billion (this will be regional for subset too)
            rgndf[!,:pct_gdp]=rgndf[!,:cost]./rgndf[!,:gdp] # Annual cost / annual GDP
            out=rgndf

        elseif sumsegs=="global"

            globdf = out |> @groupby({_.ens,_.time, _.level, _.category}) |>
                @map(merge(key(_),{cost = sum(_.cost)})) |> DataFrame

            globsoc = innerjoin(temp4,temp5,on=[:time,:ciam_country])
            globsoc[!,:gdp] = globsoc[!,:ypcc].*globsoc[!,:pop]./1e3
            globsoc[!,:ens] = fill(ensnum, size(globsoc)[1])
            globsoc = globsoc |> @groupby({_.ens,_.time}) |>
                @map(merge(key(_), {ypcc=sum(_.gdp), pop=sum(_.pop),gdp=sum(_.gdp)})) |> DataFrame
            globsoc[!,:ypcc] = (globsoc[!,:gdp].*1e9)./(globsoc[!,:pop].*1e6)

            globdf = innerjoin(globdf,globsoc,on=[:ens,:time])
            globdf[!,:segID] = fill(0., size(globdf)[1])
            globdf[!,:ciam_country] = fill("global", size(globdf)[1])
            globdf[!,:segments] = fill("global", size(globdf)[1])
            globdf=globdf[!, [:ens,:time,:ciam_country,:segments,:segID,:category,:level,:cost,:ypcc,:pop,:gdp]]
            subset[!,:gdp]=subset[!,:ypcc].*subset[!,:pop]./1e3
            globdf=[globdf;subset]

            globdf[!,:pct_gdp] = globdf[!,:cost] ./ globdf[!,:gdp]

            if rgns!=false
                rgndf = filter(row -> row[:ciam_country] in rgns,out)
                rgndf = rgndf |> @groupby({_.ens,_.time,_.ciam_country, _.level, _.category,_.ypcc,_.pop}) |> @map(merge(key(_),{cost = sum(_.cost)})) |> DataFrame
                rgndf[!,:segID] = fill(0., size(rgndf)[1])
                rgndf[!,:segments] = fill("regional", size(rgndf)[1])
                rgndf = rgndf[!,col_order]

                rgndf[!,:gdp]=rgndf[!,:ypcc].*rgndf[!,:pop]./1e3 # GDP is in $Billion (this will be regional for subset too)
                rgndf[!,:pct_gdp]=rgndf[!,:cost]./rgndf[!,:gdp] # Annual cost / annual GDP

                globdf=[globdf; rgndf]
            end

            out=globdf
        else
            out[!,:gdp]= out[!,:ypcc] .* out[!,:pop] ./1e3 # Regional gdp by segment
            out[!,:pct_gdp] = out[!,:cost] ./ out[!,:gdp] # Segment cost or subcost as % of regional gdp
        end

        if i==1
            global df = out
        else
            df = [df;out]
        end

    end

    # Remove ypcc and pop from final df
    df = df[!,[:ens,:time,:ciam_country,:segments,:segID,:category,:level,:cost,:gdp,:pct_gdp]]
    return df

end

## ---------------------------------
## 4. In Progress Helper Functions
## ---------------------------------

"""
    writelog()

IN PROGRESS - Create a writelog function to go with a wrapper for the run function
automatically produce a logfile
"""
function writelog()
    dir=joinpath(@__DIR__,"..","data","batch","logs")
    d = init()
    run = d["run_name"]
    date = Dates.now()
    cp("../data/batch/init.csv", joinpath(dir,"$(run)_$(date).csv"))
end

using CSV

## ---------------------------------
## 5. Write-Out Helper Functions
## ---------------------------------

"""
    write_init_file(run_name::String, outputdir::String, init_settings::Dict)

Write the init.csv file for a specificied `run_name` into `outputdir` using init_settings
found in `init_settings`.
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

## ---------------------------------
## 6. SLRCost Component Helper Functions
## ---------------------------------

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
