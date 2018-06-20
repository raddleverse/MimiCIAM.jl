# # Catherine Ledna
# 4/12/18
#------------------------------------------------------------------------
# ciamhelpers.jl
#------------------------------------------------------------------------
# Assorted functions to process CIAM data and run CIAM model
#------------------------------------------------------------------------

using Distributions

# Function to load CIAM parameters from CSV to dictionary
#   data_dir = relative path to data location
function load_ciam_params(data_dir)
    files = readdir(data_dir)
    filter!(i->(i!="desktop.ini" && i!=".DS_Store" && i!="xsc.csv"), files)
    params = Dict{Any, Any}(lowercase(splitext(m)[1]) => readdlm(joinpath(data_dir,m), ',' ) for m in files)
    return params

end

# Function to read in LSLR from file and filter to desired set of segments
#   Modifies input parameter dictionary
# subset - list of segments you want
# params - parameter dictionary you want to add lslr to 
function preplsl!(data_dir,lslfile,subset, params)
    lsl_params = Dict{Any, Any}("lslr" => readdlm(joinpath(data_dir,lslfile), ',' ))

    # Filter LSL according to subset segments
    p = lsl_params["lslr"]
    p_new = p[2:end,:]
    row_order = sortperm(p_new[:,1])
    p_new = p_new[row_order,1:end]

    if subset != false
        s = p_new[:,1]
        ind_s = filter_index(s,subset)
        p_new = p_new[ind_s,2:end]
    
        params["lslr"] = p_new
    else
        params["lslr"] = p_new[:, 2:end]
    end
    

    return params
end

# Function to process CIAM data from csv to usable format
#   Stores outputs in params
# rgn_order, seg_order - alphabetized lists of regions/segments used
# Specific to CIAM data; some hard-coded names and assumptions
function parse_ciam_params!(params, rgn_order, seg_order)
    key = [k for k in keys(params)]

    for i in 1:length(key)
        p = params[key[i]] # Array
        keyname = key[i]

        if keyname=="data"
            rownames = copy(p[1:1,:]) # Preserve first row names
            
            # Sort Segments alphabetically
            segnames= filter(f -> f!="NA",p[:,1])
            row_order = sortperm(segnames)
            p = p[2:end,:]
            p = p[row_order,:]

            segs = p[:,1]

            # Filter segments to desired values based on inputs
            seg_inds = filter_index(segs, seg_order)

            if length(seg_inds)>=1
                # Process all variables
                for k in 2:size(p,2)
                    varname = rownames[1,k]
                    newvars = p[1:end,k]
                    newvars = newvars[seg_inds]
                    newvars = [convert(Float64,v) for v in newvars]     
                    params[varname] = newvars
                end
                delete!(params, "data")
            else
                error("Segments in dictionary do not match supplied segments") 
            end
        elseif keyname=="globalparams"

            for k in 1:size(p,1)
                varname = p[k,1]
                newval = p[k,2]

                if (varname=="ntsteps" || varname=="adaptPers")
                    newval = convert(Int64,newval)

                elseif typeof(newval)==Int
                    newval = convert(Float64,newval)

                elseif isa(newval,AbstractString)
                    if lowercase(newval)=="true"
                        newval = true

                    elseif lowercase(newval)=="false"
                        newval = false
                    end
                end                
                params[varname] = newval

            end
            delete!(params, "globalparams")
 
        elseif size(p,2) ==2
            # Alphabetize
            p = p[sortperm(p[:,1]),:]

            # Filter regions
            r = p[:,1]
            r_inds = filter_index(r, rgn_order)
            p = p[r_inds,:]

            if p[:,1]!=rgn_order
                error("Regions in dictionary do not match supplied regions, ", keyname)               
            else
                newvals = p[:,2]
                params[keyname] = newvals
            end
        elseif size(p,2)>3
            # Country-time data matrices
            # Alphabetize
            p = p[2:end,:]  # cut off time 
            row_order = sortperm(p[:,1])
            p = p[row_order, :]

            # Filter regions and trim name column
            r = p[:,1]
            ind = filter_index(r,rgn_order)
            p = p[ind, 2:end] 

            params[keyname] = p
           
        elseif size(p,2)==1 && typeof(p)==Array{Float64,2}
            p_new = p[:,1]
            params[keyname] = p_new
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

function transpose_string_matrix(mat)
    # Utility function to manually transpose a matrix composed of non-numeric types
    nmat = vec(mat[1,:])
    for i in collect(2:size(mat,1))
        col = mat[i,:]
        nmat = hcat(nmat,col)
    end 
    return nmat   
end

# Function to process the segment-country mapping file (xsc) in CIAM
#   Reads from CSV, outputs list of dictionaries and arrays
#   Filters xsc file to desired segments/regions
function prepxsc(data_dir, xscfile, subset)

    xsc_name = replace(xscfile, ".csv","") # Strip ".csv" from file name

    # Read in csv and convert to dictionary format 
    xsc_params = Dict{Any, Any}(lowercase(splitext(xscfile)[1]) => readdlm(joinpath(data_dir, xscfile), ',' ))
    xsc_char = Dict{Any,Any}( xsc_params[xsc_name][i,1] => (xsc_params[xsc_name][i,2],xsc_params[xsc_name][i,3]) for i in 1:size(xsc_params[xsc_name],1))

    # If only a subset of segments is used, filter down to relevant segments
    if subset!=false
        filter!((k,v)-> k in subset, xsc_char)
    end

    # Create region and segment indices
    rgns = sort(unique([i[1] for i in collect(values(xsc_char))]))
    segs = sort(unique(collect(keys(xsc_char))))

    xsc_ind = Dict{Any,Any}()      # numeric seg -> (numeric rgn, greenland bool)
    xsc_segmap = Dict{Any,Any}()   # Numeric seg/rgn -> char seg/rgn
    xsc_rgnmap = Dict{Any,Any}()
   
    for i in 1:length(segs)
        r = xsc_char[segs[i]][1]   # Region character
        grn = xsc_char[segs[i]][2] # 0 = non-Greenland, 1 = greenland bool 
        r_ind = findind(r, rgns)   # Region index 
        
        new_val = (r_ind, grn)     # New tuple w/ region index instead of character
        
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
    name_ind = find(h, vec)[1]
    return name_ind

end




