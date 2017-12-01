# Helper functions from FUND (github.com/davidanthoff/fund.jl/src/helper.jl)
# Modified for CIAM by Catherine Ledna, November 10, 2017

using Distributions

function loadparametersciam(datadir=joinpath(dirname(@__FILE__), "..", "data"))
    files = readdir(datadir)
    filter!(i->(i!="desktop.ini" && i!=".DS_Store" && i!="xsc.csv"), files)
    parameters = Dict{Any, Any}(lowercase(splitext(file)[1]) => readdlm(joinpath(datadir,file), ',' ) for file in files)

    #prepparameters!(parameters)

    return parameters
end

function parse_ciam_params!(params, rgn_order, seg_order)

    # 1. String to tuple
    for i in params
        p = i[2] # Array
        keyname = i[1]

        # Special case for seg-data csv
        if keyname=="data"
            segnames= filter(f -> f!="NA",p[:,1])
            row_order = sortperm(segnames)
            segs = segnames[row_order]

            if segs == seg_order
                for k in collect(2:size(p,2))
                    varname = p[1,k]
                    newvars = p[2:size(p,1),k]
                    newvars = newvars[row_order]
                    params[varname] = newvars
                end
                delete!(testparams, "data")
            else
                error("Segments in dictionary do not match supplied segments") 
            end
 
        elseif size(p,2) ==2
            # Type 1: Country/Segment-Data Format
            # Assumption: column 1 = country column
            # Make and sort tuples; confirm region/segment order and remove region name

            p_tup = sort([ (p[j,1], p[j,2]) for j in collect(1:size(p,1))])

            if collect([j[1] for j in p_tup])!=rgn_order
                error("Regions in dictionary do not match supplied regions")                
            else
                newvals = [ j[2] for j in p_tup ]
                params[keyname] = newvals
            end

        elseif size(p,1)>3
            # Type 2: Time-Country-Data format
            # CASES: 
            #   1. Country-time matrix
            #   2. GSL case

            # Sort by country, time period
            # ASSUMPTION: Row 1 = strings of country names
            # ASSUMPTION: Column 1 = time 1-20
            if (size(p, 2)==length(rgn_order) | size(p,2)==(length(rgn_order)+1))
                
                p_new = p[:,2:size(p,2)]
                col_order = sortperm(p_new[1,:])
                p_new = p_new[2:size(p,1),col_order]
            
                params[keyname] = p_new
            end

        end

    end

end


function prepxsc(params)  
    # Returns index of regions (rgns_full), index of segments (segs_full), dictionary translating segment index
    #   to region index (xsc_ind_full), and dictionary translating segment name to region name (xsc_char)

    # Create segment-region mapping (string version)
    xsc_char = Dict{Any,Any}( split(params["xsc"][i],",")[1] =>  split(params["xsc"][i],",")[2] for i in collect(1:length(params["xsc"])))

    # Create region and segment indices
    rgns_full = unique(collect(values(xsc_char)))
    segs_full = unique(collect(keys(xsc_char)))
    
    # Map segment index to region index
    xsc_ind_full = Dict()
    for i in collect(1:length(segs_full))
        r = xsc_char[segs_full[i]]
        r_ind = findind(r, rgns_full)
        xsc_ind_full[i] = r_ind
    end
    return (xsc_ind_full, rgns_full, segs_full, xsc_char)
end

function findind(val, vec)
    # Look up index corresponding to name
    # vec - a vector of region or segment names (strings)
    # val - a string corresponding to value in 'vec'
    h(u) = u == val
    name_ind = find(h, vec)[1]
    return name_ind

end

import StatsBase.mode
function mode(d::Truncated{Gamma{Float64},Continuous})
    return mode(d.untruncated)
end

function getindexfromyear(year)
    const baseyear = 1950
    return year - baseyear + 1
end

function convertparametervalue(pv)
    if isa(pv,AbstractString)
        if startswith(pv,"~") & endswith(pv,")")
            args_start_index = search(pv,'(')
            dist_name = pv[2:args_start_index-1]
            args = split(pv[args_start_index+1:end-1], ';')
            fixedargs = filter(i->!contains(i,"="),args)
            optargs = Dict(split(i,'=')[1]=>split(i,'=')[2] for i in filter(i->contains(i,"="),args))

            if dist_name == "N"
                if length(fixedargs)!=2 error() end
                if length(optargs)>2 error() end

                basenormal = Normal(parse(Float64, fixedargs[1]),parse(Float64, fixedargs[2]))

                if length(optargs)==0
                    return basenormal
                else
                    return Truncated(basenormal,
                        haskey(optargs,"min") ? parse(Float64, optargs["min"]) : -Inf,
                        haskey(optargs,"max") ? parse(Float64, optargs["max"]) : Inf)
                end
            elseif startswith(pv, "~Gamma(")
                if length(fixedargs)!=2 error() end
                if length(optargs)>2 error() end

                basegamma = Gamma(parse(Float64, fixedargs[1]),parse(Float64, fixedargs[2]))

                if length(optargs)==0
                    return basegamma
                else
                    return Truncated(basegamma,
                        haskey(optargs,"min") ? parse(Float64, optargs["min"]) : -Inf,
                        haskey(optargs,"max") ? parse(Float64, optargs["max"]) : Inf)
                end
            elseif startswith(pv, "~Triangular(")
                triang = TriangularDist(parse(Float64, fixedargs[1]), parse(Float64, fixedargs[2]), parse(Float64, fixedargs[3]))
                return triang
            else
                error("Unknown distribution")
            end
        elseif pv=="true"
            return true
        elseif pv=="false"
            return false
        elseif endswith(pv, "y")
            return parse(Int, strip(pv,'y'))
        else
            try
                return parse(Float64, pv)
            catch e
                error(pv)
            end
        end
        return pv
    else
        return pv
    end
end

function getbestguess(p)
    if isa(p, ContinuousUnivariateDistribution)
        return mode(p)
    else
        return p
    end
end

function prepparameters!(parameters)
    for i in parameters
        p = i[2]
        column_count = size(p,2)
        if column_count == 1
            parameters[i[1]] = getbestguess(convertparametervalue(p[1,1]))
        elseif column_count == 2
            parameters[i[1]] = Float64[getbestguess(convertparametervalue(p[j,2])) for j in 1:size(p,1)]
        elseif column_count == 3
            length_index1 = length(unique(p[:,1]))
            length_index2 = length(unique(p[:,2]))
            new_p = Array(Float64,length_index1,length_index2)
            cur_1 = 1
            cur_2 = 1
            for j in 1:size(p,1)
                new_p[cur_1,cur_2] = getbestguess(convertparametervalue(p[j,3]))
                cur_2 += 1
                if cur_2 > length_index2
                    cur_2 = 1
                    cur_1 += 1
                end
            end
            parameters[i[1]] = new_p
        end
    end
end