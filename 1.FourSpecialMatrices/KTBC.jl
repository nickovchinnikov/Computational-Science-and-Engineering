module KTBC
  using Test

  import ToeplitzMatrices.Toeplitz
  import LinearAlgebra.Tridiagonal

  function CreateKTBC(dims=2)
    toeplitzArgs = [2, -1, zeros(dims - 2)...]
    K = Tridiagonal(Toeplitz(toeplitzArgs, toeplitzArgs))
    T = copy(K)
    T[1, 1] = 1
    B = copy(T)
    B[dims, dims] = 1
    C = convert(Array{Float64}, K)
    C[1, dims] = -1
    C[dims, 1] = -1
    return K, T, B, C
end

  K, T, B, C = CreateKTBC(3)

  @test K == [
    2.0 -1.0 0.0
    -1.0 2.0 -1.0
    0.0 -1.0 2.0
  ]

  @test T == [
    1.0 -1.0 0.0
    -1.0 2.0 -1.0
    0.0 -1.0 2.0
  ]

  @test B == [
    1.0 -1.0 0.0
    -1.0 2.0 -1.0
    0.0 -1.0 1.0
  ]

  @test C == [
    2.0 -1.0 -1.0
    -1.0 2.0 -1.0
    -1.0 -1.0 2.0
  ]
end
