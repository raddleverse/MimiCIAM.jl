using Mimi

cd("Desktop/ERG/Research/ciam/mimi-ciam.jl/src")
include("ciam.jl")
include("slrcost.jl")
include("ciamhelper.jl")
using Main.ciam 

run(ciam.getciam)

# Write model results to data frame ('output/results-jl')ff
write_ciam(ciam) # Write segment-level results
write_ciam(ciam,sumsegs="rgn") # Write results summed to region
write_ciam(ciam,sumsegs="global") # Write results summed to global 
