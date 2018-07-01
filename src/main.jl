include("ciamhelper.jl")

lslfile = "lsl_rcp0_p50.csv"
subset = false  # use string of subset file location relative to main.jl to run for a subset of coastal segments

m = get_ciam()  # Run ciam for SLR of rcp0p50 and all 12148 coastal segments