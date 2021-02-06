using Test

import LinearAlgebra.I
import LinearAlgebra.pinv
import LinearAlgebra.lu

include("../../../lib/KTBC.jl")
import .KTBC.CreateKTBC

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
