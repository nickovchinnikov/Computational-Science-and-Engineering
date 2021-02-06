using Test

import LinearAlgebra.I
import LinearAlgebra.det
import LinearAlgebra.pinv
import LinearAlgebra.lu

include("../../../lib/KTBC.jl")
import .KTBC.CreateKTBC

# 14
K₃, = CreateKTBC(3)

# Calculate value
function makeSingularReduceDiagonal(M, step=.1, detV=.1)
    n, m = size(M)
    if (n != m)
        throw(DimensionMismatch("M must be square matrix!"))
    end
    while true
        M = M - I(n) * step;
        det(M) <= detV && break
    end
    return M
end

@test_throws DimensionMismatch makeSingularReduceDiagonal([
  1 2 3
  1 3 4
])

singularK₃ = makeSingularReduceDiagonal(K₃, .0000001, .0000001);
singularK₃

# 1.4142134996579785 <- very close!
singularK₃[1,1]
det(singularK₃)

# Solve the equation and give the exact answer √2
M = [
  √2 -1 0
  -1 √2 -1
  0 -1 √2
]

M[1,1]

det([
  √2 -1 0
  -1 √2 -1
  0 -1 √2
])

u = [
  1
  √2
  1
]

@test round.(Int, M * u) == [
  0
  0
  0
]

# 17
# Cat't solve in case of dependent cols. Only linear regression possible!
A = [
  2 3 1
  2 3 1
  3 2 1
]
invA = pinv(A)
b = [
  1
  2
  0
]

x = invA * b

A * x

@test A * x != b

# 18

A = [
  2 3
  4 5
]
B = [
  1 2
  2 4
]

A * B

# 21
A = [
  2 0
  3 1
]

B = [
  1 2
  4 5
]

@test (A * B) * (A * B) == (A * B)^2
