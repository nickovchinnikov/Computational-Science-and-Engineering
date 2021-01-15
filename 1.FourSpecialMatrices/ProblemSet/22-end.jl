using Test
using Plots

import SparseArrays.spdiagm
import LinearAlgebra.I
import ToeplitzMatrices.Toeplitz

# 22
n = 1000
e = ones(n)

K = spdiagm(-1 => -e[1:999], 0 => e * 2, 1 => -e[1:999])

@test K[1,1] == 2.0
@test K[2,1] == -1.0
@test K[1,2] == -1.0
@test K[500,456] == 0.0
@test K[75,22] == 0.0
@test K[1000,999] == -1.0
@test K[1000,1000] == 2.0
@test K[999,1000] == -1.0

u = K \ e
plot(u)

# 23
A = rand(Float64, 4, 3)

K₂ = spdiagm(-1 => A[:,1], 0 => A[:,2], 1 => A[:,3])

# last 2 items is left out -> last(A[:,1]) and last(A[:,3])
Kₓ = convert(Array{Float64}, K₂)

@test Kₓ[5,4] == last(A[:,1]) && Kₓ[4,5] == last(A[:,3])

# 24
@test spdiagm(0 => ones(100)) == I(100)

# 26
v = [2, -1, zeros(2)..., -1]
Toeplitz(v, v)

# 27
A₀ = [
  -1 1 0 0
  0 -1 1 0
  0 0 -1 1
]

K = A₀ * A₀'

A₁ = A₀[:,2:4]

T = A₁ * A₁'

A₂ = A₀[:,2:3]

B = A₂ * A₂'
