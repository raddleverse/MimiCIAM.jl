# Catherine Ledna
# 3/7/18
#------------------------------------------------------------------------
# Test Code
#------------------------------------------------------------------------
# Run model and compare results of Julia and GAMS code; 
# Outputs csv with rank by accuracy (1-6 scale, 6 best)
# In progress for multi-segment comparisons 
#------------------------------------------------------------------------

using Mimi
include("../src/ciam.jl")
include("../src/ciamhelper.jl")
include("test.jl")
 
#------------------------------------------------------------------------
# CIAM - Single-Segment version (Philippines10615 segment)
#------------------------------------------------------------------------

data_dir_phl = joinpath("test_phil/input-data")
gamsfile = "../results-gams/test.csv" 
lslfilephl = "lsl_rcp0_p50.csv"
jlfile = "../results-jl/results.csv"
xscfile = "xsc.csv"
resultsdirphl = "test_phil/comparison"

n = run_tests(data_dir,gamsfile,jlfile,resultsdir,lslfile,["Philippines10615"],["rcp0_p50"], "ciam") # TODO: broken

# Test of Write_Results
xsc = prepxsc(data_dir, xscfile, ["Philippines10615"])
j = write_results(n[1], "rcp0_p50", ".", xsc)

# Testing w/o println
l = run_and_write_model(data_dir_phl, ".",lslfilephl,["Philippines10615"],"rcp0_p50",true,false)

#------------------------------------------------------------------------
# CIAM - Multisegment version (1000 segments)
#------------------------------------------------------------------------
data_dir = "../data/input-data"
resultsdir = "../data/results"
lslfile = "lsl_rcp85_p50.csv"
gamsfile = "../data/input-data/results/results-gams/globalCostB.csv"
jlfile = "../data/input-data/results/results-jl/results.csv"

segnames = readlines(open("../test-data-delavane/segmentnames.csv"))

# Test and Benchmark Only: No Results comparison
d = import_model_data(data_dir, lslfile, "xsc.csv", segnames)
k = run_and_write_model(data_dir, resultsdir,lslfile,segnames,"rcp85_p50", false, true)

for i in 1:20
    @time model_driver(data_dir, lslfile, segnames)
end
# Averages 1.9-2.2 seconds, 21.83 M allocations; 530.734 MiB; 22-31% gc time per run
h = model_driver(data_dir, lslfile, segnames)
xsc2 = prepxsc(data_dir, xscfile, segnames)
xsc2[4]
#------------------------------------------------------------------------
# CIAM - Full version (~12000 segments)
#------------------------------------------------------------------------

for i in 1:10
    @time model_driver(data_dir, lslfile, false)
end

# Full segments: 12.4-13.9 seconds each; 217.63 M allocatons; 4.672 GiB; 23.21% gc time)