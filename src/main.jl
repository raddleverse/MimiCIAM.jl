using Mimi
using MimiCIAM

m = MimiCIAM.get_model()
run(m)

# Write model results to data frame ('output/results-jl')
MimiCIAM.write_ciam(m) # Write segment-level results
MimiCIAM.write_ciam(m; sumsegs="rgn") # Write results summed to region
MimiCIAM.write_ciam(m; sumsegs="global") # Write results summed to global 
