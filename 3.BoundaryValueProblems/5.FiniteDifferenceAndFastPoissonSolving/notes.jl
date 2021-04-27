using Test
using LinearAlgebra

import Pkg
Pkg.add("FFTW")
using FFTW

Fₘ = fft(I(4), 1)

imag(Fₘ)

function CreateLM(N=2)
  L = 2 * ones(1, N) - 2 * cos.(1:N)' * π / (N + 1)
  LM = ones(N, 1) * L + L' * ones(1, N)
end

CreateLM(2)

# Packing 2D Vector into Matrices 
N = 4
L = 2 * ones(1, N) - 2 * cos.(1:N)' * π / (N + 1)
LM = ones(N, 1) * L + L' * ones(1, N)
