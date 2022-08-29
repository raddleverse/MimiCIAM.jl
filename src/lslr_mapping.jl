"""
Replaces the followingi n the event you want tight coupling

    preplsl!(lslfile,subset, params,segnames)

Read in LSLR from file and filter to desired set of segments, note that this modifies
the input parameter dictionary `params`.  The arguments are as fullows:

- lslfile - name of lslr file to use; location relative to data/input-data directory
- subset - list of segments you want
- params - parameter dictionary you want to add lslr to
- segnames - names of segments

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

using Mimi

@defcomp lsl_mapping begin

    # --- Indices ---
    segments = Index()

    # --------------------
    # Model Parameters
    # --------------------

#    lslr = Parameter(index = [time, segments], unit = "m")                # Local sea level rise (m)
    # slr_gsic
    # slr_gis
    # slr_ais
    # slr_te
    # slr_lws

    # --------------------
    # Model Variables
    # --------------------

    lslr = Variable(index = [time, segments], unit = "m")                # Local sea level rise (m)

    # --------------------
    # Model Equations
    # --------------------

    function run_timestep(p, v, d, t)

        # This is a workaround for a type instability that should be fixed in Mimi.jl
        d_segments = d.segments::Vector{Int}

        if is_first(t)
            for m in d_segments
                v.lslr[t,m] = 0.0
            end
        else
            for m in d_segments
                v.lslr[t,m] = 0.0
            end
        end

    end
end
