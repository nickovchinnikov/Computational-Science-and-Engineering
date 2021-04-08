using LinearAlgebra
using Test

# 1
# a
A = [
   3 -1 -1  0 -1
  -1  3  0 -1 -1
  -1  0  3 -1 -1
   0 -1 -1  3 -1
  -1 -1 -1 -1  4
]

AᵀA = A' * A

# b
@test tr(AᵀA) == 68

@test rank(AᵀA) == 4

values, = eigen(AᵀA)

@test round(sum(values)) == 68

@test round(reduce(*, values)) == 0

# c
A₁ = [
   1 -1 0 0 0
  -1  1 0 0 0
   0  0 0 0 0
   0  0 0 0 0
   0  0 0 0 0
]

A₁ᵀA₁ = A₁' * A₁

@test A₁ᵀA₁ == [
  2  -2  0  0  0
 -2   2  0  0  0
  0   0  0  0  0
  0   0  0  0  0
  0   0  0  0  0
]

#d
b = ones(8, 1)

@test_throws SingularException inv(A)
