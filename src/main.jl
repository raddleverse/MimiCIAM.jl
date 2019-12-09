using Mimi

include("ciam.jl")
include("slrcost.jl")
include("ciamhelper.jl")
using Main.MimiCIAM

m = MimiCIAM.get_model()
run(m[1])

# Write model results to data frame ('output/results-jl')
write_ciam(m) # Write segment-level results
write_ciam(m,sumsegs="rgn") # Write results summed to region
write_ciam(m,sumsegs="global") # Write results summed to global 
