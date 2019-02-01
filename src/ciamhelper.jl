# # Catherine Ledna
# 4/12/18
#------------------------------------------------------------------------
# ciamhelpers.jl
#------------------------------------------------------------------------
# Assorted functions to process CIAM data and run CIAM model
#------------------------------------------------------------------------
include("slrcost.jl")

using Mimi
using Distributions
using DelimitedFiles
using DataFrames
using Query
using CSV


# Function to load CIAM parameters from CSV to dictionary
function load_ciam_params()
    data_dir = "../data/input"
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
    data_dir = "../data/lslr"
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

    data_dir = "../data/input"
    xscfile = "xsc.csv"
    xsc_name = replace(xscfile, r".csv" => s"") # Strip ".csv" from file name

    # Read in csv and convert to dictionary format 
    xsc_params = Dict{Any, Any}(lowercase(splitext(xscfile)[1]) => CSV.read(joinpath(data_dir, xscfile)))
    xsc_char = Dict{Any,Any}( xsc_params[xsc_name][i,1] => (xsc_params[xsc_name][i,2],xsc_params[xsc_name][i,3], xsc_params[xsc_name][i,4]) for i in 1:size(xsc_params[xsc_name],1))

    # If only a subset of segments is used, filter down to relevant segments
    if subset!=false
        filter!((k,v)-> k in subset, xsc_char)
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
   
    return (xsc_ind, rgns, segs, xsc_char, xsc_rgnmap, xsc_segmap)

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
    dir="../data/subsets"
    if subset!=false
        subs=readlines(joinpath(dir,subset))
        return subs
    else
        return false
    end
end

function init()
    dir="../data/meta"
    header = DelimitedFiles.readdlm(joinpath(dir,"init.txt"),',')
    lsl=header[1]
    subset=header[2]
    return(lsl,subset)

end

# Wrapper for importing model data.
# lslfile - filename for lsl (string)
# subset - filename with names of segments to use (string) or false (bool) to run all segments
function import_model_data()
    # Extract lslfile and subset from data
    (lslfile,sub)=init()
    subset=load_subset(sub)

    # Process main and lsl params
    params = load_ciam_params()
     
    # Process XSC (segment-country mapping dictionary)
    xsc = prepxsc(subset)

    # Process params using xsc and format lsl file
    parse_ciam_params!(params, xsc[2], xsc[3])
    preplsl!(lslfile, subset, params,xsc[3])

    return(params, xsc)

end

function load_meta()
    metadir = "../data/meta"

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
# m - an instance of the model
# RCP - string for RCP we're using; todo make automatic from lsl file
# xsc - segment-region dictionaries
# outputdir, outfile - where to write results to, relative to test folder
# sumsegs - whether to sum across segments
function write_ciam(m, xsc; rcp="na", tag="untagged", sumsegs=false)
    outputdir = "../output/results-jl"
    meta_output = load_meta()

    outfile = "results_$(rcp)_$(tag).csv"
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
        d = m[:ciam, Symbol(v)]
            
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
                val = d[:,i] 
                v_arr = func(val)
    
                if level=="protect"
                    for j in 1:4
                        p = protectdict[j]
                        lev_arr = repeat([p], outer = n)
    
                        outarr= [rcp_arr lev_arr s_arr cost_arr t v_arr]
                                # Write results to csv 
                 
                        open(joinpath(outputdir,outfile),"a") do g 
                            writedlm(g, outarr,',')
                        end
    
                    end
                elseif level=="retreat"
                    for j in 1:5
                        q = retreatdict[j]
                        lev_arr = repeat([q], outer = n)
    
                        outarr= [rcp_arr lev_arr s_arr cost_arr t v_arr]  
                
                        open(joinpath(outputdir,outfile),"a") do g 
                            writedlm(g, outarr,',')
                        end
                    end
    
                else
                    lev_arr = repeat([level], outer = n)
                    outarr= [rcp_arr lev_arr s_arr cost_arr t v_arr]
    
                    open(joinpath(outputdir,outfile),"a") do g 
                        writedlm(g, outarr,',')
                    end
    
                end
    
            elseif typeof(d)==Array{Float64,3}
                for j in 1:size(d,3)
                    val = d[:,i,j] 
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
                        writedlm(g, outarr,',')
                    end
    
                end
            end
        end
    end
end



