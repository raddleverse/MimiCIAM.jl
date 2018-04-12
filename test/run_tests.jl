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
jlfile = "../results-jl/results.csv"
resultsdir = "test_phil/comparison"

n = run_tests(data_dir,gamsfile,jlfile,resultsdir,["rcp0_p50"], "ciam")

#------------------------------------------------------------------------
# CIAM - Multisegment version (1000 segments)
#------------------------------------------------------------------------


