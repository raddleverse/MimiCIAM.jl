include("ciamhelper.jl")

lslfile = "lsl_rcp0_p50.csv"
subset = false  # use string of subset file location relative to main.jl to run for a subset of coastal segments

m = get_ciam(lslfile, subset)  # Run ciam for SLR of rcp0p50 and all 12148 coastal segments

write_ciam(m[1], m[2], rcp="rcp0p50", tag="allsegs", sumsegs=false) # Write ciam results to output/results-jl directory
