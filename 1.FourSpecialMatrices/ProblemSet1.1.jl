using Test

import LinearAlgebra.diagm
import LinearAlgebra.triu
import LinearAlgebra.I
import LinearAlgebra.det

include("./KTBC.jl")
import .KTBC.CreateKTBC

# 3
U_5 = I(5) - diagm(1 => ones(4))
S_5 = triu(ones(5, 5))

U_5 * S_5

@test S_5 == inv(U_5)

# 4
S_4 = triu(ones(4, 4))

S_4 * S_4'

# 5
KTBC_2 = CreateKTBC()
K_2 = convert(Array{Int}, KTBC_2[1])

det(K_2)
inv(K_2)

@test inv(K_2) * det(K_2) == [
  2.0 1.0
  1.0 2.0
]

KTBC_3 = CreateKTBC(3)
K_3 = convert(Array{Int}, KTBC_3[1])

det(K_3)
inv(K_3)

@test round.(Int, inv(K_3) * det(K_3)) == [
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
