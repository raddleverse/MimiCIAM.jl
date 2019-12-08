
function buildciam(m::Model)
    add_comp!(m, slrcost)
    return m

end

function initciam(xsc, params, m::Model, t::Int=20)

    discountrate = 0.04#parse(Float64,d["discountrate"])
    
    # Dynamically find indices corresponding to USA and CAN and manually set time steps 
    rgn_ind_canada = [k for (k,v) in xsc[5] if v=="CAN"][1]
    rgn_ind_usa = [k for (k,v) in xsc[5] if v=="USA"][1]

    # Add component: slrcost and set some parameters manually 
    set_param!(m, :slrcost, :xsc, xsc[1])
    set_param!(m, :slrcost, :rgn_ind_canada, rgn_ind_canada)
    set_param!(m, :slrcost, :rgn_ind_usa, rgn_ind_usa)
    set_param!(m, :slrcost, :discountrate, discountrate)
    set_param!(m, :slrcost, :ntsteps, t)

    # Shorten some time-dependent parameters to correspond to the correct number of timesteps 
    for k in keys(params)
        if size(params[k])!=() && length(size(params[k]))==2 && size(params[k])[1]>t && k!="surgeexposure"
            params[k]=params[k][1:t,:]
        end
    end

    # Set the rest of the parameters 
    set_leftover_params!(m, params)

end

function get_model(t::Int=20)
    initparams=init()
    modelparams = import_model_data(initparams["lslr"][1],initparams["subset"][1])
    params = modelparams[1]
    xsc = modelparams[2]

    run_name = initparams["run_name"]

    
    m=Model()

    set_dimension!(m, :time, t)
    set_dimension!(m, :adaptPers, length(params["at"]))  
    set_dimension!(m, :regions, xsc[2])
    set_dimension!(m, :segments, xsc[3])

    buildciam(m)
    initciam(xsc, params, m, t)

    return (m,xsc,initparams) 

end

