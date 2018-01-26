# Compare results with Delavane's model

# 1. Import data 
# 2. 

using Mimi
using ExcelReaders
include("../src/ciam.jl")
include("../src/helper.jl")

# Temp: Convert Delavane data to usable params
# lslr = readxlsheet("../test-data-delavane/lslr.xlsx", "Sheet2")
# data = readxlsheet("../test-data-delavane/datasub1.xlsx", "Sheet2")
# writecsv("../test-data-delavane/csvs/lslr.csv", lslr)
# writecsv("../test-data-delavane/csvs/data.csv", data)
# gdp = readxlsheet("../test-data-delavane/CIAM gdp pop.xlsx", "ym")
# gdp = gdp[3:4082,:]
# writecsv("../test-data-delavane/csvs/ym.csv", gdp)


# Data preparation
datadir = joinpath("..", "test-data-delavane/csvs")
parameters = loadparametersciam(datadir)

file="xsc.csv"
xsc_params = Dict{Any, Any}(lowercase(splitext(file)[1]) => readdlm(joinpath(datadir,file), '\r' ))
rgn_seg_params = prepxsc(xsc_params)
xsc_ind_full = rgn_seg_params[1]
allrgns = sort(rgn_seg_params[2])

allsegs = sort(rgn_seg_params[3])
xsc_char_full = rgn_seg_params[4]
parse_ciam_params!(parameters, allrgns, allsegs)

# Post-processing for pop and GDP 
pop_table = parameters["pop"]
rgns = pop_table[1,:][2:161]
rgn_order = sortperm(rgns)
rgn_ordered = rgns[rgn_order]

# Want: array of indices from rgn_ordered corresp to positions of members of allrgns
function findind(val, vec)
    # Look up index corresponding to name
    # vec - a vector of region or segment names (strings)
    # val - a string corresponding to value in 'vec'
    h(u) = u == val
    name_ind = find(h, vec)[1]
    return name_ind

end
outind = [0]
for i in allrgns
    ind = findind(i, rgn_ordered)
    if outind==[0]
        outind = ind
    else
        outind = vcat(outind,ind)
    end
end

pop_ordered = pop_table[:,2:161]
pop_ordered = pop_ordered[:, rgn_order]
pop_ordered = pop_ordered[:, outind]
pop_ordered = pop_ordered[2:21,:]

gdp_table = parameters["ypcc"][:, 2:161]
rgnorder = sortperm(gdp_table[1,:])
gdp_ordered = gdp_table[:, rgnorder]
gdp_ordered = gdp_ordered[:, outind][2:21,:]

ym = parameters["ym"][:,2:205]
rgnorder = sortperm(ym[1,:])
ym_ord = ym[:, rgnorder]

rgns2 = ym[1,:]
rgnord2 = sortperm(rgns2)
rgnsord2 = rgns2[rgnord2]

outind2 = [0]
for i in allrgns
    ind = findind(i, rgnsord2)
    if outind2==[0]
        outind2 = ind
    else
        outind2 = vcat(outind2,ind)
    end
end


ym_ord = ym_ord[:, outind2]
ym_ord = ym_ord[2:21,:]

parameters["ypcc"] = gdp_ordered
parameters["pop"] = pop_ordered
#parameters["ym"] = 

# -----------------------------------------------------
# Run model 
atpers = [1, 5, 10]
m = Model()
setindex(m, :time, 14)
setindex(m, :adaptPers, length(atpers))
setindex(m, :regions, allrgns)
setindex(m, :segments, allsegs)

addcomponent(m, ciam)

# Time params
setparameter(m, :ciam, :ntsteps, 14)       
setparameter(m, :ciam, :tstep, 1.)
setparameter(m, :ciam, :at, atpers) 

# Adaptation Parameters
adapt = [1,10,100,1000,10000]
setparameter(m, :ciam, :adaptOptions, adapt)
setparameter(m, :ciam, :fixed, true)

# To do: NA handling
# For popdens: if not given by segment, then use country scaling 
#   To do: Implement this in CIAM 
# To do: Account for other NAs and CIAM implementation

# Today's discussion: 
# Add methods for when some strings are NA
# Send 


setleftoverparameters(m, parameters)

