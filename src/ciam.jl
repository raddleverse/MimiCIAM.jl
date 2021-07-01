using Mimi

"""
    initciam(xsc, params, initparams, m::Model; fixed::Bool = false, t::Int = 20, 
            noRetreat::Bool = false, allowMaintain::Bool = true, popinput::Int = 0)

Initialize a CIAM model `m` with the given arguments.
"""
function initciam(xsc, params, initparams, m::Model; fixed::Bool = false, t::Int = 20, 
                noRetreat::Bool = false, allowMaintain::Bool = true, popinput::Int = 0)

    discountrate = 0.04

    if initparams["lslr"][1] !="lsl_rcp0_p50.csv"
        rcp = parse(Int64,replace(replace(initparams["lslr"][1],r"^[^l]*lsl_rcp"=>s""),r"_.*"=>s""))
        pctl = parse(Int64,replace(replace(initparams["lslr"][1], r"^[^l]*lsl_rcp[0-9][0-9]_p"=>s""),r".csv"=>s""))
    else
        rcp = 0
        pctl = 50
    end

    if initparams["ssp"][1] == false
        ssp=0
    else
        ssp = parse(Int64,replace(replace(initparams["ssp"][1],r"^[^l]*SSP"=>s""),r"_.*"=>s""))
    end

    # Dynamically find indices corresponding to USA and CAN and manually set time steps
    # If the lengths are 0, then assume those segments are not used. Note that
    # if including Greenland, need Canada too as a reference for land appreciation
    
    rgn_ind_canada = [k for (k,v) in xsc[4] if v=="CAN"]
    rgn_ind_canada = (length(rgn_ind_canada) > 0) ? rgn_ind_canada[1] : 0

    rgn_ind_usa = [k for (k,v) in xsc[4] if v=="USA"]
    rgn_ind_usa = (length(rgn_ind_usa) > 0) ? rgn_ind_usa[1] : 0

    # Add component: slrcost and set some parameters manually
    segID = segStr_to_segID(xsc[3])
    
    update_param!(m, :slrcost, :segID, segID)
    update_param!(m, :slrcost, :rcp, rcp)
    update_param!(m, :slrcost, :percentile, pctl)
    update_param!(m, :slrcost, :ssp, ssp)
    update_param!(m, :slrcost, :xsc, xsc[1])
    update_param!(m, :slrcost, :rgn_ind_canada, rgn_ind_canada)
    update_param!(m, :slrcost, :rgn_ind_usa, rgn_ind_usa)
    update_param!(m, :slrcost, :discountrate, discountrate)
    update_param!(m, :slrcost, :ntsteps, t)
    update_param!(m, :slrcost, :fixed, fixed)
    update_param!(m, :slrcost, :noRetreat, noRetreat)
    update_param!(m, :slrcost, :allowMaintain, allowMaintain)
    update_param!(m, :slrcost, :popinput, popinput)

    # Shorten some time-dependent parameters to correspond to the correct number of timesteps
    for k in keys(params)
        if ndims(params[k]) == 0 && ndims(params[k]) == 2 && size(params[k])[1] > t && k != "surgeexposure"
            params[k] = params[k][1:t,:]
        end
    end

    # Set the rest of the parameters - to use update_leftover_params! we need To
    # specific the component as well as the parameter name
    params_dict = Dict()
    for (k,v) in params
        params_dict[(:slrcost, k)] = v
    end
    update_leftover_params!(m, params_dict)

end

"""
    get_model(;initfile::Union{String, Nothing} = nothing, fixed::Bool=false, 
                t::Int = 20, noRetreat::Bool = false, allowMaintain::Bool = true, 
                popinput::Int = 0)
Return a initialized and built CIAM model with the given arguments.
"""
function get_model(;initfile::Union{String, Nothing} = nothing, fixed::Bool=false, 
                    t::Int = 20, noRetreat::Bool = false, allowMaintain::Bool = true, 
                    popinput::Int = 0)

    initparams  = init(; f = initfile)
    params, xsc = import_model_data(initparams["lslr"][1], initparams["subset"][1], initparams["ssp"][1], initparams["ssp_simplified"][1])

    m = Model()

    set_dimension!(m, :time, t)
    set_dimension!(m, :adaptPers, length(params["at"]))
    set_dimension!(m, :regions, xsc[2])
    set_dimension!(m, :segments, xsc[3])

    add_comp!(m, slrcost)

    initciam(xsc, params, initparams, m; fixed = fixed, t = t, noRetreat = noRetreat, allowMaintain = allowMaintain, popinput = popinput)

    return m

end

# TODO - Code to run a batch instance of the model
# function update_model(updateparams,updatevalues)
# end
