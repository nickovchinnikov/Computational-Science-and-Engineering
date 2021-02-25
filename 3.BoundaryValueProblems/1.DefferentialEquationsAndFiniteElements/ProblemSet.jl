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
