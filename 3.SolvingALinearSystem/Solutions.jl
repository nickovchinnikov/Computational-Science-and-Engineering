using Test
using Plots

include("../lib/KTBC.jl")
using .KTBC: CreateKTBC

using LinearAlgebra: I, SymTridiagonal, lu, det, Diagonal, diagm

# 1
K₄ = CreateKTBC(4)[1]

Lₐ = [
  1 0 0 0
  -1 / 2 1 0 0
  0 -2 / 3 1 0
  0 0 -3 / 4 1
]

Uₐ = [
  2 -1 0 0
  0 3 / 2 -1 0
  0 0 4 / 3 -1
  0 0 0 5 / 4
]

L, U = lu(K₄)

@test K₄ == L * U
@test K₄ == L * Diagonal(U) * L'
@test det(K₄) == 5

# 3
K₅ = CreateKTBC(4)[1]

@test round.(det(K₅) * inv(K₅)) == [
 4.0  3.0  2.0  1.0
 3.0  6.0  4.0  2.0
 2.0  4.0  6.0  3.0
 1.0  2.0  3.0  4.0
]

L, U = lu(K₅)
invL = inv(L)

@test K₅ == L * Diagonal(U) * L'

function checkInvL(i::Int, j::Int)
    return round(invL[i, j], digits=2) == round(j / i, digits=2)
end

@test checkInvL(3, 1)
@test checkInvL(3, 2)
@test checkInvL(4, 2)
@test round(invL[4, 2], digits=1) == 1 / 2

# 4
d = (2:5)./(1:4)

L = I(4) - diagm(-1 => (1:3)./(2:4))

L * Diagonal(d) * L'
