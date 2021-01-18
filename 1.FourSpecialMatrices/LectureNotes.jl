using Test

using ToeplitzMatrices: Toeplitz
using SparseArrays: sparse
using LinearAlgebra: triu, diag, diagm, Diagonal, Bidiagonal, Tridiagonal, det, I

include("../lib/KTBC.jl")
using .KTBC: CreateKTBC

# Difference matrix
K₄ = Toeplitz([2, -1, 0, 0], [2, -1, 0, 0])
K₁₀ = Toeplitz([2, -1, zeros(10)...], [2, -1, zeros(10)...])

# T - Top boundary condition
RectangularT = Toeplitz([1, -1, 0, 0], [1, 0, 0])

@test RectangularT == [
  1.0   0.0   0.0
 -1.0   1.0   0.0
  0.0  -1.0   1.0
  0.0   0.0  -1.0
]

# Circulant C, periodic
C₄ = Toeplitz([2, -1, 0, -1], [2, -1, 0, -1])

@test C₄ == [
  2.0  -1.0   0.0  -1.0
 -1.0   2.0  -1.0   0.0
  0.0  -1.0   2.0  -1.0
 -1.0   0.0  -1.0   2.0
]

nullspaceVector = C₄ * ones(4)

@test nullspaceVector == [0.0; 0.0; 0.0; 0.0]

# Sparse version of K₄
sparse(K₄)

invK₄ = round.(inv(convert(Array{Float64}, K₄)); digits=1)

@test invK₄ == [
 0.8  0.6  0.4  0.2
 0.6  1.2  0.8  0.4
 0.4  0.8  1.2  0.6
 0.2  0.4  0.6  0.8
]

U = triu(ones(3, 3))

@test U == [
 1.0  1.0  1.0
 0.0  1.0  1.0
 0.0  0.0  1.0
]

function createDiagonal(value=1, dim=4, isMatrix=true)
    DiagonalMatrix = zeros(Int, dim, dim) + I * value
    return isMatrix ? DiagonalMatrix : diag(DiagonalMatrix)
end

# Not memory efficient method
function createKWithMemoryProblem(dim=4)
    DiagonalMatrix = createDiagonal(2, dim)
    E = diagm(1 => -ones(dim - 1))
    K = DiagonalMatrix + E + E'
end

@test createKWithMemoryProblem() == [
  2.0  -1.0   0.0   0.0
 -1.0   2.0  -1.0   0.0
  0.0  -1.0   2.0  -1.0
  0.0   0.0  -1.0   2.0
]

# Very efficient
function createK(dim=4)
    ULDiagonals = ones(Int, dim - 1)
    Tridiagonal(ULDiagonals, createDiagonal(2, dim, false), ULDiagonals)
end

convert(Array{Int}, createK())

@test convert(Array{Int}, createK()) == [
 2  1  0  0
 1  2  1  0
 0  1  2  1
 0  0  1  2
]

# Convert to memory-efficient
Tridiagonal(K₁₀)

K, T, B, C = CreateKTBC(8)

@test [det(K), det(T), det(B), round(det(C))] == [9.0; 1.0; 0.0; 0.0]

# Worked examples
e = ones(8)
u = [1:8...]

@test B * e == [
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
]

f = B * u

@test f' * e == 0.0

# First way to make fixed-free
H = convert(Array{Int}, Toeplitz([2, -1, 0], [2, -1, 0]))
H[3, 3] = 1

@test H == [
  2  -1   0
 -1   2  -1
  0  -1   1
]

# Second way
J = reverse(convert(Array{Int}, I(8)), dims=2)
T = convert(Array{Int}, T)
H = J * T * J

# Inverse ways
invH = round.(Int, inv(H))
invHUsingT = round.(Int, J * inv(T) * J)

@test invH == invHUsingT
