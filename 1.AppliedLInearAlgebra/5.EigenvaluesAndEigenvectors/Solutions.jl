using Test
using LinearAlgebra

include("../../lib/KTBC.jl")
using .KTBC: CreateKTBC

# 1
K₂ = CreateKTBC(2)[1]

Q = (1 / √2) * [
  1 1
  1 -1
]

Λ = [
  1 0
  0 3
]

@test K₂ == round.(Q * Λ * Q')

# 3, 4
K₅ = convert(Array{Int}, CreateKTBC(5)[1])
eK₅ = round.(eigvals(K₅), digits=3)

@test eK₅ == round.([
  2 - √3
  2 - 1
  2 - 0
  2 + 1
  2 + √3
], digits=3)

@test eK₅ == round.(2 * ones(5, 1) - 2 * cos.((1:5) * π / 6), digits=3)[1:5]

Q = round.(eigvecs(K₅), digits=3)
DST_C = round.(Q * diagm(0 => [1,-1,-1,-1,-1]), digits=3)

JK = (1:5) * (1:5)'

DST = sin.(JK * π / 6)/√3

@test DST_C == round.(DST, digits=3)
@test round.(DST', digits=3) == round.(DST^-1, digits=3)

# 6
T₆ = convert(Array{Int}, CreateKTBC(6)[2])
eT₆ = round.(eigvals(T₆), digits=3)

@test eT₆ == round.(2 * ones(6, 1) - 2 * cos.((.5:5.5) * π / 6.5), digits=3)[1:6]

# 7
C₄ = convert(Array{Int}, CreateKTBC(4)[4])
Λ, Q = eigen(C₄)

F = [
  1 1 1 1
  1 im im^2 im^3
  1 im^2 im^4 im^6
  1 im^3 im^6 im^9
]

A = Q \ F

# 16
A₁ = [
  .6 .4
  .4 .6
]

eigvals(A₁)

@test !all(1 .> eigvals(A₁))

# (A₂)ᵏ -> 0
A₂ = [
  .6 .9
  .1 .6
]

eigvals(A₂)

@test all(1 .> eigvals(A₂))

# 19
B = [
  3 1
  0 2
]

evalsB, evecB = eigen(B)

@test B^10 == evecB * (diagm(0 => evalsB)^10) * evecB^-1

