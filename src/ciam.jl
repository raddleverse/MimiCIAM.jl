
function buildciam(m::Model)
    add_comp!(m, slrcost)
    return m

end

function initciam(xsc, params, initparams, m::Model, fixed::Bool=false,t::Int=20)

    discountrate = 0.04#parse(Float64,d["discountrate"])
    if initparams["lslr"][1] !="lsl_rcp0_p50.csv"
        rcp = parse(Int64,replace(replace(initparams["lslr"][1],r"^[^l]*lsl_rcp"=>s""),r"_.*"=>s""))
        pctl = parse(Int64,replace(replace(initparams["lslr"][1], r"^[^l]*lsl_rcp[0-9][0-9]_p"=>s""),r".csv"=>s""))
    else
        rcp=0
        pctl=50
    end
    
    if initparams["ssp"][1]==false
        ssp=0
    else
        ssp = parse(Int64,replace(replace(initparams["ssp"][1],r"^[^l]*SSP"=>s""),r"_.*"=>s""))
    end

    # Dynamically find indices corresponding to USA and CAN and manually set time steps 
    rgn_ind_canada = [k for (k,v) in xsc[4] if v=="CAN"][1]
    rgn_ind_usa = [k for (k,v) in xsc[4] if v=="USA"][1]

    # Add component: slrcost and set some parameters manually 
    segID = segStr_to_segID(xsc[3])
    set_param!(m, :slrcost, :segID, segID)
    set_param!(m, :slrcost, :rcp, rcp)
    set_param!(m, :slrcost, :percentile, pctl)
    set_param!(m, :slrcost, :ssp, ssp)
    set_param!(m, :slrcost, :xsc, xsc[1])
    set_param!(m, :slrcost, :rgn_ind_canada, rgn_ind_canada)
    set_param!(m, :slrcost, :rgn_ind_usa, rgn_ind_usa)
    set_param!(m, :slrcost, :discountrate, discountrate)
    set_param!(m, :slrcost, :ntsteps, t)
    set_param!(m, :slrcost, :fixed, fixed)

    # Shorten some time-dependent parameters to correspond to the correct number of timesteps 
    for k in keys(params)
        if size(params[k])!=() && length(size(params[k]))==2 && size(params[k])[1]>t && k!="surgeexposure"
            params[k]=params[k][1:t,:]
        end
    end

    # Set the rest of the parameters 
    set_leftover_params!(m, params)

end

function get_model(;initfile=nothing,fixed::Bool=false,t::Int=20)
    initparams= init(initfile)
    modelparams = import_model_data(initparams["lslr"][1],initparams["subset"][1],initparams["ssp"][1])
    params = modelparams[1]
    xsc = modelparams[2]

    m=Model()

    set_dimension!(m, :time, t)
    set_dimension!(m, :adaptPers, length(params["at"]))  
    set_dimension!(m, :regions, xsc[2])
    set_dimension!(m, :segments, xsc[3])

    buildciam(m)

    initciam(xsc, params, initparams, m, fixed, t)

    return m 

end

# Code to run a batch instance of the model 
function update_model(updateparams,updatevalues)

end