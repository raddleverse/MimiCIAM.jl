using Mimi

cd("Desktop/ERG/MastersProject/pjt/ciam/src/")
include("ciam.jl")
include("ciamhelper.jl")
using Main.ciam 

run(ciam.m)

# Write model results to data frame ('output/results-jl')
write_ciam(ciam) # Write segment-level results
write_ciam(ciam,sumsegs="rgn") # Write results summed to region
write_ciam(ciam,sumsegs="global") # Write results summed to global 
