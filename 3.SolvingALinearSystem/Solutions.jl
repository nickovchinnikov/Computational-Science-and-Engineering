using Test
using Plots

include("../lib/KTBC.jl")
using .KTBC: CreateKTBC

using LinearAlgebra: I, SymTridiagonal, lu, cholesky, det, eigen, Diagonal, diag, diagm, PosDefException

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
d = (2:5) ./ (1:4)

L = I(4) - diagm(-1 => (1:3) ./ (2:4))

L * Diagonal(d) * L'

# 9
K₃, T₃, B₃, C₃ = CreateKTBC(3)

A = cholesky(convert(Array{Float64}, K₃))
@test round.(A.L * A.U) == K₃

A = cholesky(convert(Array{Float64}, T₃))
@test round.(A.L * A.U) == T₃

@test_throws PosDefException cholesky(convert(Array{Float64}, B₃))

Bₑ = B₃ + eps() * I(3)
A = cholesky(convert(Array{Float64}, Bₑ))
@test round.(A.L * A.U) == round.(Bₑ)

# 10
eye₄ = ones(4, 4)
eigen(eye₄)

@test det(eye₄) == 0.0

L = [
  1
  1
  1
  1
]

@test L * L' == eye₄

# 11
K = ones(4, 4) + I(4) / 100
L, U = lu(K)
eigen(K)
inv(K)

# 12 
# https://rosettacode.org/wiki/Pascal_matrix_generation#Julia
function pascal(n::Int)
    return [binomial(j + i, i) for i in 0:n, j in 0:n]
end

# https://en.wikipedia.org/wiki/Pascal_matrix
function pascalL(n::Int)
    return exp(diagm(-1 => 1:(n - 1)))
end

K = pascal(4)
L, U = lu(K)

@test det(K) == abs(reduce(*, diag(U)))

L = pascalL(5)

@test K == round.(L * L')

# 13
Fib = [
  1 1
  1 0
]

v = [
  1
  0
]

@test Fib^6 * v == [
  13
  8
]
