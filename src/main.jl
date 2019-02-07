using Mimi

cd("Desktop/ERG/MastersProject/pjt/ciam/src/")
include("ciam.jl")
include("slrcost.jl")
include("ciamhelper.jl")
using Main.ciam 

run(ciam.m)

# Write model results to data frame ('output/results-jl')
write_ciam(ciam.m, ciam.modelparams[2]) # Write segment-level results
write_ciam(ciam.m,ciam.modelparams[2],sumsegs="rgn") # Write results summed to region
write_ciam(ciam.m,ciam.modelparams[2],tag="fullv5",sumsegs="global") # Write results summed to global 
