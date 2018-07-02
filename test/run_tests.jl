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

gamsfilephl = "../output/results-gams/Phl_test_comp.csv" 
lslfilephl = "lsl_rcp0_p50.csv"

segnames_phl = ["Philippines10615"]
rcp_phl = "rcp0_p50"

test1seg = run_tests(gamsfilephl, lslfilephl, segnames_phl, rcp_phl,"Phl")

#------------------------------------------------------------------------
# CIAM - Multisegment version (1000 segments)
#------------------------------------------------------------------------
lslfile_1000_ref = "lsl_rcp0_p50.csv"
gamsfile_1000 = "../output/results-gams/globalCostBcomp.csv"
rcp_1000_ref = "rcp0_p50"

seg1 = readlines(open("../data/subsets/sub1.csv"))
seg2 = readlines(open("../data/subsets/sub2names.csv"))
seg3 = readlines(open("../data/subsets/sub3names.csv"))
seg4 = readlines(open("../data/subsets/sub4names.csv"))
seg5 = readlines(open("../data/subsets/sub5names.csv"))
seg6 = readlines(open("../data/subsets/sub6names.csv"))
seg7 = readlines(open("../data/subsets/sub7names.csv"))
seg8 = readlines(open("../data/subsets/sub8names.csv"))
seg9 = readlines(open("../data/subsets/sub9names.csv"))
seg10 = readlines(open("../data/subsets/sub10names.csv"))


# tsub1 = run_tests(gamsfile_1000, lslfile_1000_ref, seg1, rcp_1000_ref,"sub1")
# tsub2 = run_tests(gamsfile_1000, lslfile_1000_ref, seg2, rcp_1000_ref,"sub2")
# tsub3 = run_tests(gamsfile_1000, lslfile_1000_ref, seg3, rcp_1000_ref,"sub3")
# tsub4 = run_tests(gamsfile_1000, lslfile_1000_ref, seg4, rcp_1000_ref,"sub4")
# tsub5 = run_tests(gamsfile_1000, lslfile_1000_ref, seg5, rcp_1000_ref,"sub5")
# tsub6 = run_tests(gamsfile_1000, lslfile_1000_ref, seg6, rcp_1000_ref,"sub6")
# tsub7 = run_tests(gamsfile_1000, lslfile_1000_ref, seg7, rcp_1000_ref,"sub7")
# tsub8 = run_tests(gamsfile_1000, lslfile_1000_ref, seg8, rcp_1000_ref,"sub8")
# tsub9 = run_tests(gamsfile_1000, lslfile_1000_ref, seg9, rcp_1000_ref,"sub9")
# tsub10 = run_tests(gamsfile_1000, lslfile_1000_ref, seg10, rcp_1000_ref,"sub10")

# Averages 1.9-2.2 seconds, 21.83 M allocations; 530.734 MiB; 22-31% gc time per runt

#------------------------------------------------------------------------
# CIAM - Full version (~12000 segments)
#------------------------------------------------------------------------
lslfile_full = "lsl_rcp0_p50.csv"
gamsfile_full = "../output/results-gams/total_mockup.csv"  # Don't have comparison data 
rcp_full = "rcp0_p50"

testAllsegs = run_tests(gamsfile_full,lslfile_full, false, rcp_full,"full", sum=true)




