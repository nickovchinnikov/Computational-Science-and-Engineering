using SymPy
using Test
using LinearAlgebra
using Plots
using Calculus

# using SymbolicUtils

include("../lib/KTBC.jl")
using .KTBC: CreateKTBC

using Pkg
Pkg.add("Pluto")
Pluto.run()

# 1 Symbolic path
@vars u₁ u₂ x₁ y₁

u⁺ = [u₁ u₂]
u = [
  u₁
  u₂
]

T₂ = convert(Array{Int}, CreateKTBC(2)[2])

u⁺T₂u = u⁺ * T₂ * u

result₁ = [u₁ * (u₁ - u₂) + u₂ * (-u₁ + 2 * u₂)]

@test u⁺T₂u[1] - result₁[1] == 0

# 3
A = [
  1 -1 0
  0 1 -1
  -1 0 1
]

plot3d(A)

C = (A') * A

plot!(C)

@test_throws PosDefException cholesky(C)

@test A[1:3] + A[4:6] + A[7:9] == [0; 0; 0]

# or
u = ones(Int, 3, 1)
z = zeros(Int, 3, 1)

@test A * u == z
