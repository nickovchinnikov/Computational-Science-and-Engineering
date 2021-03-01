using Test

using LinearAlgebra

include("../../lib/KTBC.jl")
import .KTBC.CreateKTBC

# 10
KTBC₃ = CreateKTBC(3)
K₃ = convert(Array{Int}, KTBC₃[1])

res = 4 * K₃

invRes = inv(res)

u(x) = x - x^2

@test u.([
  1/4
  1/2
  3/4
]) == round.(invRes * (1 / 2) * [
  1
  1
  1
], digits=4)

# 11
K = [
    12 -3
  -3   12
]
inv(K) * 1/3 * [
  1
  2
]

f(x) = - (x^3)/6 + x/6

[
  f(1/3)
  f(2/3)
]
