using Test
using LinearAlgebra

include("./KTBC.jl")
using .KTBC: CreateKTBC

function CreateK2D(dimSquared=2)
  IdentityFixedDims = I(dimSquared)
  K, T, B, C = CreateKTBC(dimSquared)
  return kron(IdentityFixedDims, K) + kron(K, IdentityFixedDims)
end

IMatrixTestDims = I(2)

@test CreateK2D() == [
   K + 2 * IMatrixTestDims   -IMatrixTestDims
  -IMatrixTestDims            K + 2 * IMatrixTestDims
]

export CreateK2D
