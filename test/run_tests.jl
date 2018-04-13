# Catherine Ledna
# 3/7/18
#------------------------------------------------------------------------
# Test Code
#------------------------------------------------------------------------
# Run model and compare results of Julia and GAMS code; 
# Outputs csv with rank by accuracy (1-6 scale, 6 best)
#------------------------------------------------------------------------

using Mimi
include("../src/ciam.jl")
include("../src/ciamhelper.jl")
include("test.jl")
 
#------------------------------------------------------------------------
# CIAM - Single-Segment version (Philippines10615 segment)
#------------------------------------------------------------------------

data_dir = joinpath("test_phil/input-data")
gamsfile = "../results-gams/test.csv" 
lslfile = "lsl_rcp0_p50.csv"
jlfile = "../results-jl/results.csv"
xscfile = "xsc.csv"
resultsdir = "test_phil/comparison"

n = run_tests(data_dir,gamsfile,jlfile,resultsdir,lslfile,["Philippines10615"],["rcp0_p50"], "ciam")

# Test of Write_Results
xsc = prepxsc(data_dir, xscfile, ["Philippines10615"])
j = write_results(n[1], "rcp0_p50", ".", xsc)

#------------------------------------------------------------------------
# CIAM - Multisegment version (1000 segments)
#------------------------------------------------------------------------
data_dir = "../data/input-data"
resultsdir = "../data/results"
lslfile = "lsl_rcp85_p50.csv"
gamsfile = "../data/input-data/results/results-gams/globalCostB.csv"
jlfile = "../data/input-data/results/results-jl/results.csv"
