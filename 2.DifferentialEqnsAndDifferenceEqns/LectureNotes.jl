module Lecture2

using Test

include("../lib/KTBC.jl")
import .KTBC.CreateKTBC

include("../lib/SpecialVectors.jl")
using .SpecialVectors: constants, linear, squares, rampAtK, sines, cosines, exponentials

include("../lib/DifferencesMatrix.jl")
using .DifferencesMatrix: backwardDiff, forwardDiff

dims = 4

K, T, B, C = CreateKTBC(dims)

Kₐ = convert(Array{Int}, K)

cV = constants(dims)

@test (K * cV)[2:3] == [
  0.0
  0.0
]

lV = linear(dims)

@test (K * lV)[2:3] == [
  0.0
  0.0
]

sV = squares(dims)

@test (-K * sV)[2:3] == [
 2.0
 2.0
]

rV = rampAtK(2, dims)

@test (-K * rV)[1:3] == [
  0.0
  1.0
  0.0
]

t = π / 4
multiplier = round(2 * cos(t) - 2, digits=2)

sinV = round.(sines(t, dims), digits=2)
@test round.((-K * sinV)[1:3], digits=1) == round.((multiplier * sinV)[1:3], digits=1)

cosV = round.(cosines(t, dims), digits=2)
@test round.((-K * cosV)[2:3], digits=1) == round.((multiplier * cosV)[2:3], digits=1)

expV = round.(exponentials(t, dims), digits=2)
@test round.((-K * expV)[2:3], digits=1) == round.((multiplier * expV)[2:3], digits=1)

# Matlab experiment
function mathlabExperiment(n::Int)
    h = 1 / (n + 1)
    u = cos.(π * [1:n...] * h / 2)
    c = (π / 2)^2
    f = c * u

    K, T, B, C = CreateKTBC(n)

    U = (h^2) * T \ f
    e = 1 - U[1]
    g = [c / 2; f]

    T = backwardDiff(n + 1) * forwardDiff(n + 1)

    V = (h^2) * T \ g
    E = 1 - V[1]

    return f, e, T, V, E
end

f, e, T, V, E = mathlabExperiment(3)

end
