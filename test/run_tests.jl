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
include("../src/helper.jl")
include("test.jl")
 
###------
# CIAM - Simplified
###------

data_dir = joinpath("test_phil/input-data")
paramfiles =  ["data.csv", "globalparams.csv","pop.csv","ypcc.csv","ypc_usa.csv"]
lsldata = "philinputlsl_reshape.csv"
gamsfile = "../results-gams/test.csv" 
jlfile = "../results-jl/results.csv"
resultsdir = "test_phil/comparison"

n = run_tests(data_dir,paramfiles,gamsfile,jlfile,resultsdir,["rcp0_p50"], "ciam")

m = run_tests(data_dir,paramfiles,gamsfile,jlfile,resultsdir,["rcp0_p50"])


d = import_comparison_data(resultsdir, "comparison.csv")


