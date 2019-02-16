module ciam

using Mimi

export m

include("slrcost.jl")
include("ciamhelper.jl")

(run_name,lsl,subset) = init()
modelparams = import_model_data(lsl,subset)
params = modelparams[1]
xsc = modelparams[2]

rgn_ind_canada = [k for (k,v) in xsc[5] if v=="CAN"][1]

m = Model()
set_dimension!(m, :time, 20)
set_dimension!(m, :adaptPers, 5)  # Todo figure out way not to hardcode these
set_dimension!(m, :regions, xsc[2])
set_dimension!(m, :segments, xsc[3])

add_comp!(m, slrcost)
set_param!(m, :slrcost, :xsc, xsc[1])
set_param!(m, :slrcost, :rgn_ind_canada, rgn_ind_canada)
set_leftover_params!(m, params)

end