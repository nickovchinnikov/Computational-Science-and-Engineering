using Test

import LinearAlgebra.diagm
import LinearAlgebra.triu
import LinearAlgebra.I
import LinearAlgebra.det

include("../KTBC.jl")
import .KTBC.CreateKTBC

function createU(dim=3)
    return I(dim) - diagm(1 => ones(dim - 1))
end

function createS(dim=3)
    return triu(ones(dim, dim))
end

# 3
U₅ = createU(5)
S₅ = createS(5)

U₅ * S₅

@test S₅ == inv(U₅)

# 4
S₄ = triu(ones(4, 4))

S₄ * S₄'

# 5
KTBC₂ = CreateKTBC()
K₂ = convert(Array{Int}, KTBC₂[1])

det(K₂)
inv(K₂)

@test inv(K₂) * det(K₂) == [
  2.0 1.0
  1.0 2.0
]

KTBC₃ = CreateKTBC(3)
K₃ = convert(Array{Int}, KTBC₃[1])

det(K₃)
inv(K₃)

@test round.(Int, inv(K₃) * det(K₃)) == [
  3 2 1
  2 4 2
  1 2 3
]

# 6
function klcounter(n=4)
    j = 1
    result = zeros(n, n)
    while j <= n
        i = j
        while i <= n
            result[i, j] = (n - i + 1) * j
            i += 1
        end
        j += 1
    end
    return result
end

@test klcounter() == [
  4.0  0.0  0.0  0.0
  3.0  6.0  0.0  0.0
  2.0  4.0  6.0  0.0
  1.0  2.0  3.0  4.0
]

# 7
KTBC₃ = CreateKTBC(3)

T₃ = convert(Array{Int}, KTBC₃[2])
K₃ - T₃

u = [
  1
  0
  0
]
vᵀ = transpose(u)

@test u * vᵀ == K₃ - T₃

invT₃ = inv(T₃)
invK₃ = inv(K₃)

@test round.(Int, 4(invT₃ - invK₃)) == [
 9  6  3
 6  4  2
 3  2  1
]

u = [
  3
  2
  1
]
vᵀ = transpose(u)

@test 1 / 4 * u * vᵀ == round.(invT₃ - invK₃; digits=2)
