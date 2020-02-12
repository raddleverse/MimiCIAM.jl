# # Catherine Ledna
# 4/12/18
#------------------------------------------------------------------------
# ciamhelpers.jl
#------------------------------------------------------------------------
# Assorted functions to process CIAM data and run CIAM model
#------------------------------------------------------------------------



# Function to load CIAM parameters from CSV to dictionary
function load_ciam_params()
    data_dir = joinpath(@__DIR__, "..","data","input")
    files = readdir(data_dir) 
    filter!(i->(i!="desktop.ini" && i!=".DS_Store" && i!="xsc.csv"), files)

    params = Dict{Any, Any}(lowercase(splitext(m)[1]) => CSV.read(joinpath(data_dir,m)) |> DataFrame for m in files)

    return params

end

# Function to read in LSLR from file and filter to desired set of segments
#   Modifies input parameter dictionary
# lslfile - string; name of lslr file to use; location relative to data/input-data directory
# subset - list of segments you want
# params - parameter dictionary you want to add lslr to 
function preplsl!(lslfile,subset, params,segnames)
    data_dir = joinpath(@__DIR__,"..","data","lslr")
    lsl_params = CSV.read(joinpath(data_dir,lslfile)) |> DataFrame

    # Filter according to subset segments
    if subset != false
        col_names = [i for i in names(lsl_params) if string(i) in subset]
        lsl_params = lsl_params[col_names]
    end

    # Chomp off unrelated rows and sort alphabetically (do this regardless of whether there's a subset)
    col_names = [i  for i in names(lsl_params) if string(i) in segnames]
    col_names = sort(col_names)
    lsl_params = lsl_params[col_names]
    
    params["lslr"] = convert(Array{Float64,2},lsl_params)

    return params
end

function prepssp!(ssp,params,rgnnames)
    data_dir = joinpath(@__DIR__,"..","data","ssp")
    if ssp==false # Do nothing, base ssp data already loaded 
        return params
    else
        pop = CSV.read(joinpath(data_dir,string("pop_",ssp,".csv")))
        ypc = CSV.read(joinpath(data_dir, string("ypcc_",ssp,".csv")))

        col_names = [i  for i in names(pop) if string(i) in rgnnames]
        col_names = sort(col_names)
        pop = pop[col_names]
        ypc = ypc[col_names]
        
        params["pop"] = Array{Float64,2}(pop)
        params["ypcc"]= Array{Float64,2}(ypc)

    end
    return params
end

# Function to process CIAM data from csv to usable format
#   Stores outputs in params
# rgn_order, seg_order - alphabetized lists of regions/segments used
# Specific to CIAM data; some hard-coded names and assumptions
function parse_ciam_params!(params, rgn_order, seg_order)
    key = [k for k in keys(params)]

    for k in key
        p = params[k] # Data frame

        if k=="data"
            colnames = filter(f -> string(f)!="NA",names(p)) # Preserve column names 
            
            # Filter segments to subset
            segs = p[1]
            seg_inds = filter_index(segs, seg_order)
            p = p[seg_inds,:]
            # Sort alphabetically
            seg_alpha = sortperm(p[1])
            p = p[seg_alpha,:]

            if length(seg_inds)>=1
                # Process all variables
                for k in 2:(length(colnames)+1)
                    varname = string(colnames[k-1])
                    newvars = p[1:end,k]
                    newvars = [convert(Float64,v) for v in newvars]     
                    params[varname] = newvars
                end
                delete!(params, "data")
            else
                error("Segments in dictionary do not match supplied segments") 
            end
        elseif k=="globalparams"

            for k in 1:length(p[1])
                varname = p[1][k]
                newval = p[2][k]

                if (varname=="ntsteps" || varname=="adaptPers")
                    newval = parse(Int64,newval)

                elseif lowercase(newval)=="true"
                    newval = true
                elseif lowercase(newval)=="false"
                    newval = false
                else
                    newval=parse(Float64,newval)
                end                
                params[varname] = newval

            end
            delete!(params, "globalparams")
        elseif k=="surgeexposure"
            p=@from i in p begin
                @where i.segments in seg_order
                @select i
                @collect DataFrame
            end
            # Sort alphabetically
            sort!(p, :segments)
            params["surgeexposure"] = convert(Array{Float64,2},p[:,2:6])
        elseif k=="refa"
            p = @from i in p begin
                @where i.segments in seg_order
                @select i
                @collect DataFrame
            end
            sort!(p,:segments)
            params["refA"]= convert(Array{Float64},p[:value])
            delete!(params,"refa")
        elseif size(p,2) ==2
            # Filter regions
            r_inds = filter_index(p[1], rgn_order)
            p = p[r_inds,:]
            # Alphabetize
            p = p[sortperm(p[1]),:]

            if p[1]!=rgn_order
                error("Regions in dictionary do not match supplied regions, ", k)               
            else
                newvals = p[2]
                # Coerce to Array{Float64,1}
                params[k] = Array{Float64,1}(newvals)
            end
        elseif size(p,2)>3
            # Time-country data matrices
            # Alphabetize
            col_names = [i for i in names(p) if string(i) in rgn_order]
            p= p[sort(col_names)]

            params[k] = Array{Float64,2}(p)
           
        elseif size(p,2)==1 
            params[k] = Array{Float64,1}(p[1])

        end  

    end

end
 
# Filters a vector (v1) by a second vector (v2), returns
#   indices of contained elements
function filter_index(v1, v2)
    out = []
    for i in 1:length(v1)
        if v1[i] in v2
            push!(out, i)
        end
    end
    return(out)
end

# Function to process the segment-country mapping file (xsc) in CIAM
#   Reads from CSV, outputs list of dictionaries and arrays
#   Filters xsc file to desired segments/regions
function prepxsc(subset)

    data_dir = joinpath(@__DIR__,"..","data","input")
    xscfile = "xsc.csv"
    xsc_name = replace(xscfile, r".csv" => s"") # Strip ".csv" from file name

    # Read in csv and convert to dictionary format 
    xsc_params = Dict{Any, Any}(lowercase(splitext(xscfile)[1]) => CSV.read(joinpath(data_dir, xscfile)))
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

# Function to look up index corresponding to name
# vec - a vector of region or segment names (strings)
# val - a string corresponding to value in 'vec'
function findind(val, vec)
    h(u) = u == val
    name_ind = findall(h, vec)[1]
    return name_ind

end

function load_subset(subset=false)
    dir=joinpath(@__DIR__,"..","data","subsets")
    if subset!=false
        subs=readlines(joinpath(dir,subset))
        return subs
    else
        return false
    end
end

function init(f=nothing)
    if f==nothing
        f=joinpath(@__DIR__,"..","data","batch","init.txt")
    end
    varnames=CSV.read(open(f),header=true)
    vardict = Dict{Any,Any}( String(i) => varnames[i] for i in names(varnames))
    return(vardict)
end

# In progress create a writelog function to go with a wrapper for the run function 
# automatically produce a logfile 
function writelog()
    dir=joinpath(@__DIR__,"..","data","batch","logs")
    d = init()
    run = d["run_name"]
    date = Dates.now()
    cp("../data/batch/init.txt",joinpath(dir,"$(run)_$(date).txt"))

end

# Wrapper for importing model data.
# lslfile - filename for lsl (string)
# subset - filename with names of segments to use (string) or false (bool) to run all segments
function import_model_data(lslfile,sub,ssp) 

    subset=load_subset(sub)

    # Process main and lsl params
    params = load_ciam_params()
     
    # Process XSC (segment-country mapping dictionary)
    xsc = prepxsc(subset)

    # Process params using xsc and format lsl file
    parse_ciam_params!(params, xsc[2], xsc[3])
    preplsl!(lslfile, subset, params,xsc[3])
    prepssp!(ssp,params,xsc[2])

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

# Function to write out model results to CSV file
# m - output from get_model() 
# sumsegs - whether to sum across all segments, to region level, or no sums
# varnames: if not false, write the passed variable names; if false get defaults from file
# To do: possibly modify to work with DataVoyager()
function write_ciam(m; runname="base", sumsegs="seg", varnames=false,tag="")
    outputdir = joinpath(@__DIR__,"..","output")
    meta_output = load_meta()
    rcp = m[:slrcost,:rcp]
    pctl = m[:slrcost,:percentile]
    rcp_str = "$(rcp)p$(pctl)"
    
    model = m
    xsc = load_xsc()
    segmap = load_segmap()
    segRgnDict = Dict{Any,Any}( xsc[:seg][i] => xsc[:rgn][i] for i in 1:size(xsc,1))

    if varnames==false
        varnames = [k for k in keys(meta_output[2])] # to do change
    end

    vargroup1 = [] # 2D vars
    vargroup2 = [] # vars greater than 2D

    for v in varnames
        if length(size(model[:slrcost,Symbol(v)]))>2
            push!(vargroup2, Symbol(v))
        else
            push!(vargroup1, Symbol(v))
        end

    end
 
    # Assign 2D variables to dataframe
    # 2 cases: 1. adapt pers is first; 2. adapt pers is second 
    common_order = [:time,:regions,:segments,:level]
    for i in 1:length(vargroup1)
        temp = getdataframe(model,:slrcost => vargroup1[i])

        missing_names = [j for j in common_order if !(j in names(temp))]
        if length(missing_names)>=1 
            temp[missing_names]=Missing
        end
        if :regions in missing_names && !(:segments in missing_names)
            temp = temp |> @map(merge(_,{regions=segRgnDict[_.segments]})) |> DataFrame
        end
        
        temp[:variable]= fill(String(vargroup1[i]),nrow(temp))
        rename!(temp,vargroup1[i]=>:value)
        temp = temp[[:time,:regions,:segments,:level, :variable, :value]]
            
        if i==1
            global df = temp
        else 
            df = [df;temp]
        end
    end

    # Assign 3D variables to second data frame and join
    ntime = model[:slrcost,:ntsteps]
    segID = model[:slrcost,:segID]
    colnames = Symbol.(segID_to_seg(Int64.(segID),segmap))

    for j in 1:length(vargroup2)
        ndim1 = size(model[:slrcost,vargroup2[j]])[3]

        for k in 1:ndim1

            temp = DataFrame(model[:slrcost,vargroup2[j]][:,:,k])

            names!(temp,colnames)
            temp[:time] = 1:ntime       
            
            if String(vargroup2[j])=="Construct" || occursin("Protect",String(vargroup2[j]))
                dim1 = k+1
                adapt=model[:slrcost,:adaptoptions][dim1]
            else
                dim1=k
                adapt=model[:slrcost,:adaptoptions][dim1]
            end

            temp[:level] = fill(adapt, ntime)
            temp = stack(temp,colnames)
            rename!(temp,:variable => :segments)
            temp[:segments] = [String(i) for i in temp[:segments]]

            temp[:variable]= fill(String(vargroup2[j]),nrow(temp))
            
            temp = temp |> @map(merge(_,{regions=segRgnDict[_.segments]})) |> DataFrame
            temp = temp[[:time,:regions,:segments,:level,:variable,:value]]


            if j==1 && k==1
                global df2 = temp
            else
                df2 = [df2;temp]
            end

        end

    end

    # Sum to either region-level, global-level, or leave as seg-level  
    outdf = [df;df2]
    outfile = "$(runname)_$(sumsegs)_$(rcp_str)_$(tag).csv"

    if sumsegs=="rgn"
        rgndf = outdf |> @groupby({_.time,_.regions, _.level, _.variable}) |> @map(merge(key(_),{value = sum(_.value)})) |> DataFrame
        CSV.write(joinpath(outputdir,outfile),rgndf)  
    elseif sumsegs=="seg"
        CSV.write(joinpath(outputdir,outfile),outdf)  
    elseif sumsegs=="global"
        globdf = outdf |> @groupby({_.time, _.level, _.variable}) |> 
            @map(merge(key(_),{value = sum(_.value)})) |> DataFrame
        CSV.write(joinpath(outputdir,outfile),globdf) 
    end  

    
end

# Function to streamline writing results for optimal adaptation costs 
function write_optimal_costs(m;runname="base")
    # Output: Data Frame with segment,region,time,level,option, suboption
    #   E.g. 'OptimalFixedProtect', 'Construct'
    # Should output 2 CSVs: 1 with just the 3 main categories, 2nd with 
    #   detailed subcategories
    outputdir = joinpath(@__DIR__,"..","output")
    rcp = m[:slrcost,:rcp]
    pctl = m[:slrcost,:percentile]
    rcp_str = "$(rcp)p$(pctl)"
    
    model = m
    xsc = load_xsc()
    segRgnDict = Dict{Any,Any}( xsc[:seg][i] => xsc[:rgn][i] for i in 1:size(xsc,1))

    # 1. Create aggregate adaptation decision DF 
    temp1 = getdataframe(model, :slrcost => :OptimalFixedCost)
    temp1 = temp1 |> @map(merge(_,{regions=segRgnDict[_.segments]})) |> DataFrame

    temp2 = getdataframe(model, :slrcost => :OptimalFixedLevel)
    temp3 = getdataframe(model, :slrcost => :OptimalFixedOption)

    # Join dataframes and reorganize
    out = join(temp1,temp2, on=[:time,:segments])
    out = join(out,temp3, on=[:time,:segments])
    
    # Replace OptimalFixedOption numeric value with string
    lookup = Dict{Any,Any}(-2.0=> "RetreatCost", -1.0=> "ProtectCost",-3.0=>"NoAdaptCost")
    out = out |> @map(merge(_,{variable=lookup[_.OptimalFixedOption]})) |> DataFrame
    rename!(out, Dict(:OptimalFixedLevel => :level))
    out = out[[:time,:regions,:segments,:variable,:level,:OptimalFixedCost]]

    # Write to file 
    outfile = "$(runname)_seg_$(rcp_str)_optimal.csv"
    CSV.write(joinpath(outputdir,outfile),out)  

    # Write Sub-Costs 
    vars = [:OptimalFixedStormCapital, :OptimalFixedStormPop, :OptimalFixedConstruct,
            :OptimalFixedFlood, :OptimalFixedRelocate, :OptimalFixedWetland]

    
    for i in 1:length(vars)
        temp = getdataframe(model, :slrcost => vars[i])
        temp = temp |> @map(merge(_,{regions=segRgnDict[_.segments]})) |> DataFrame

        temp[:variable]= fill(String(vars[i]),nrow(temp))

        temp2 = getdataframe(model, :slrcost => :OptimalFixedLevel)
        temp3 = getdataframe(model, :slrcost => :OptimalFixedOption)

        # Join dataframes and reorganize
        out = join(temp,temp2, on=[:time,:segments])
        out = join(out,temp3, on=[:time,:segments])

        # Replace OptimalFixedOption numeric value with string
        lookup = Dict{Any,Any}(-2.0=> "RetreatCost", -1.0=> "ProtectCost",-3.0=>"NoAdaptCost")
        out = out |> @map(merge(_,{AdaptCategory=lookup[_.OptimalFixedOption]})) |> DataFrame
        rename!(out, Dict(:OptimalFixedLevel => :level))
        rename!(out,vars[i]=>:value)
        out = out[[:time,:regions,:segments,:AdaptCategory,:variable,:level,:value]]

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

function write_optimal_protect_retreat(m; runname="base")
    outputdir = joinpath(@__DIR__,"..","output")
    rcp = m[:slrcost,:rcp]
    pctl = m[:slrcost,:percentile]
    rcp_str = "$(rcp)p$(pctl)"

    model = m
    xsc = load_xsc()
    segRgnDict = Dict{Any,Any}( xsc[:seg][i] => xsc[:rgn][i] for i in 1:size(xsc,1))

    # 1. Create aggregate adaptation decision DF 
    pl = getdataframe(model, :slrcost => :OptimalProtectLevel)
    rl = getdataframe(model, :slrcost => :OptimalRetreatLevel)
    out = join(pl,rl, on=[:time,:segments]) 

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

### Basic Functions for Segment-Region Lookup
function load_segmap()
    segmap = CSV.read(joinpath(@__DIR__, "..", "data","meta", "segIDmap.csv")) |> DataFrame
    segmap[:segID] = [Int64(i) for i in segmap[:segID]]
    return(segmap)
end

function load_xsc()
    xsc = CSV.read(joinpath(@__DIR__,"..","data","input","xsc.csv")) |> DataFrame
    return(xsc)
end

# Look up string name of segment from segID
# segID = int or array of ints
# segmap = output of load_rgnmap or load_segmap() (DataFrame)
# Returns only first result for each ID entry as an array 
function segID_to_seg(segID, segmap)
    seg = [String(segmap[:seg][segmap.segID.==i][1]) for i in segID]
    return(seg)
end

# Get segID from string name of segment
# Segstr must be an array of strings 
function segStr_to_segID(segstr)
    ids = [parse(Int64,replace(i, r"[^0-9]"=> "")) for i in segstr] 
    return(ids)
end