# Catherine Ledna
# 3/7/18
#------------------------------------------------------------------------
# Test Code
#------------------------------------------------------------------------
# Run model and compare results of Julia and GAMS code
#------------------------------------------------------------------------

using Mimi
include("../src/ciam.jl")
include("../src/ciamhelper.jl")
include("test.jl")
 
#------------------------------------------------------------------------
# CIAM - Single-Segment version (Philippines10615 segment)
#------------------------------------------------------------------------

data_dir_phl = "../data/input-data"
gamsfilephl = "../results/results-gams/Phl_test_comp.csv" 
lslfilephl = "lsl_rcp0_p50.csv"
resultsdirphl = "../data/results/results-jl"

segnames_phl = ["Philippines10615"]
rcp_phl = "rcp0_p50"

test1seg = run_tests(data_dir_phl, gamsfilephl, resultsdirphl, lslfilephl, segnames_phl, rcp_phl)

#------------------------------------------------------------------------
# CIAM - Multisegment version (1000 segments)
#------------------------------------------------------------------------
data_dir_1000 = "../data/input-data"
resultsdir_1000 = "../data/results/results-jl"
lslfile_1000 = "lsl_rcp85_p50.csv"
gamsfile_1000 = "../results/results-gams/globalCostBcomp.csv"

segnames_1000 = readlines(open("../test-data-delavane/segmentnames.csv"))
rcp_1000 = "rcp85_p50"

test1000segs = run_tests(data_dir_1000, gamsfile_1000, resultsdir_1000, lslfile_1000, segnames_1000, rcp_1000)

# Averages 1.9-2.2 seconds, 21.83 M allocations; 530.734 MiB; 22-31% gc time per run

#------------------------------------------------------------------------
# CIAM - Full version (~12000 segments)
#------------------------------------------------------------------------
data_dir_full = "../data/input-data"
resultsdir_full = "../data/results/results-jl"
lslfile_full = "lsl_rcp85_p50.csv"
gamsfile_full = "../results/results-gams/total_mockup.csv"  # Don't have comparison data 
rcp_full = "rcp85_p50"

testAllsegs = run_tests(data_dir_full, gamsfile_full, resultsdir_full, lslfile_full, false, rcp_full)


#for i in 1:10
#    @time model_driver(data_dir, lslfile, false) # TODO re-write this function 
#end




