using Test

include("../../lib/DifferencesMatrix.jl")
using .DifferencesMatrix: backwardDiff, forwardDiff, centeredDiff

# 8
dim = 5
forwardDiff(dim) - backwardDiff(dim)

backwardDiff(dim) * forwardDiff(dim)

backwardDiff(dim) * centeredDiff(dim)

# 10
forwardDiff(dim) * backwardDiff(dim)
