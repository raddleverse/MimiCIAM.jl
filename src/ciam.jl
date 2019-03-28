module ciam

using Mimi

export m

include("slrcost.jl")
include("ciamhelper.jl")

# Load initialization parameters: run name, lsl source, segments to run, SSP, discountrate 
d = init()
run_name = d["run_name"]
modelparams = import_model_data(d["lslr"],d["subset"],d["ssp"])
params = modelparams[1]
xsc = modelparams[2]
discountrate = parse(Float64,d["discountrate"])

# Dynamically find indices corresponding to USA and CAN and manually set time steps 
rgn_ind_canada = [k for (k,v) in xsc[5] if v=="CAN"][1]
rgn_ind_usa = [k for (k,v) in xsc[5] if v=="USA"][1]
t=10 # Running 10 periods (thru 2100)

# Set model dimensions 
m = Model()
set_dimension!(m, :time, t)
set_dimension!(m, :adaptPers, 1)  # Todo figure out way not to hardcode these
set_dimension!(m, :regions, xsc[2])
set_dimension!(m, :segments, xsc[3])

# Add component: slrcost and set some parameters manually 
add_comp!(m, slrcost)
set_param!(m, :slrcost, :xsc, xsc[1])
set_param!(m, :slrcost, :rgn_ind_canada, rgn_ind_canada)
set_param!(m, :slrcost, :rgn_ind_usa, rgn_ind_usa)
set_param!(m, :slrcost, :discountrate, discountrate)

# Shorten some time-dependent parameters to correspond to the correct number of timesteps 
for k in keys(params)
    if size(params[k])!=() && length(size(params[k]))==2 && size(params[k])[1]>t
        params[k]=params[k][1:t,:]
    end
end

# Set the rest of the parameters 
set_leftover_params!(m, params)

end