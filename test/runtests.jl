using Mimi
using MimiCIAM
using Test

include("utils.jl")

run_baseline_comparisons = true # slow so only include it when have the time etc.

include("test_unit_tests.jl")
run_baseline_comparisons && include("test_baseline_comparisons.jl")
