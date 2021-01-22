using Test
using Plots

include("../../lib/DifferencesMatrix.jl")
using .DifferencesMatrix: backwardDiff, forwardDiff, centeredDiff

# 8
dim = 5
forwardDiff(dim) - backwardDiff(dim)

backwardDiff(dim) * forwardDiff(dim)

backwardDiff(dim) * centeredDiff(dim)

# 10
forwardDiff(dim) * backwardDiff(dim)

# 13
function discretSolution(n::Int)
    h = 1 / (n + 1)
    T = backwardDiff(n) * (-forwardDiff(n))
    u = h^2 * inv(T) * ones(n)
    return u
end

trueSolutions(x) = (1 / 2) * (1 - x^2)

steps = 0:1/8:1

uₜ = trueSolutions.([steps...])[2:8]
plot(uₜ)

u₇ = discretSolution(7)
plot!(u₇)

e = uₜ - u₇
plot(e)
