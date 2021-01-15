using Test

import LinearAlgebra.diagm
import LinearAlgebra.triu
import LinearAlgebra.I
import LinearAlgebra.det
import LinearAlgebra.pinv
import LinearAlgebra.lu

include("./KTBC.jl")
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

# 8
K₅, T₅ = CreateKTBC(5)

6 * (inv(T₅) - inv(K₅))

6 * inv(K₅)

# 9
K₄, T₄, B₄, C₄  = CreateKTBC(4)
f = [
  1
  -1
  1
  -1
]
e = [
  1
  1
  1
  1
]

u = pinv(C₄) * f

function roundWithLosts(value, lost=10^(-5))
    round.(Int, lost * value)
end

# left division operator: \ - action its just like inv(a) * b 
@test roundWithLosts(C₄ \ e) == roundWithLosts(inv(C₄) * e)
@test roundWithLosts(C₄ \ f) == roundWithLosts(inv(C₄) * f)

# 10
H₃ = [
  2 -1 0
  -1 2 -1
  0 -1 1
]

invH₃ = round.(Int, inv(H₃))

I₃ = convert(Array{Int}, I(3))

J₃ = reverse(I₃; dims=2)

@test J₃ * T₃ * J₃ == H₃

@test invH₃ == J₃ * round.(Int, inv(T₃)) * J₃

U₃ = createU(3)
S₃ = createS(3)

@test U₃ * U₃' == H₃

@test S₃' * S₃ == invH₃

# 11

U₃ = createS(3)

SE = [
 0.0  0.0  1.0
 0.0  1.0  1.0
 1.0  1.0  1.0
]

NW = [
 1.0  1.0  1.0
 1.0  1.0  0.0
 1.0  0.0  0.0
]

L = [
 1.0  0.0  0.0
 1.0  1.0  0.0
 1.0  1.0  1.0
]

@test J₃ * U₃ == SE

@test U₃ * J₃ == NW

@test J₃ * U₃ * J₃ == L

@test SE * NW == [
 1.0  0.0  0.0
 2.0  1.0  0.0
 3.0  2.0  1.0
]

# 12, 13
L, U = lu(C₄)

@test round.(Int, U * 6) == [
  12  -6   0  -6
  0   9  -6  -3
  0   0   8  -8
  0   0   0   0
]

U = 6 \ round.(Int, U * 6)

# Only non-zero rows
U = U[vec(mapslices(row -> any(row .!= 0), U', dims=1)), :]
Uᵀ = U'

UUᵀ = U * Uᵀ

UUᵀ = 36 \ round.(Int, 36 * U * Uᵀ)
U⁺ = Uᵀ * inv(UUᵀ)

L = C₄ * U⁺

# Proof
@test C₄ == round.(Int, L * U)

# Or 
Lₐ = C₄ * pinv(U)
@test C₄ == round.(Int, Lₐ * U)

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
