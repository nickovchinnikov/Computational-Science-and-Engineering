using Test
using Plots

include("../../lib/DifferencesMatrix.jl")
using .DifferencesMatrix: backwardDiff, forwardDiff, centeredDiff

include("../../lib/SpecialVectors.jl")
using .SpecialVectors: squares

include("../../lib/KTBC.jl")
using .KTBC: CreateKTBC

# 8
dim = 5
forwardDiff(dim) - backwardDiff(dim)

backwardDiff(dim) * forwardDiff(dim)

backwardDiff(dim) * centeredDiff(dim)

# 10
forwardDiff(dim) * backwardDiff(dim)

# 13
function discretSolution(n::Int)
    h = 1 / (n + 1)
    T = backwardDiff(n) * (-forwardDiff(n))
    u = h^2 * inv(T) * ones(n)
    return u
end

n = 7
analogSolution(x) = (1 / 2) * (1 - x^2)
steps = 0:1 / (n + 1):1

uₜ = analogSolution.([steps...])[2:(n + 1)]
plot(uₜ)

u₇ = discretSolution(n)
plot!(u₇)

e = uₜ - u₇
plot(e)

# 14
function cosAgainstDiscretCos(n::Int)
    x = 0:10^(-2):1

    h = 1 / (n + 1)
    xᵢ = 0:h:1
    ω = 4 * π

    u = cos.(ω * x)
    uᵢ = cos.(ω * h * (1:n + 2))

    lineWidth = 2
    plot(x, u, label="analog cos", lw=lineWidth, title="Cos/DiscretCos Oscilations with N=$n")
    plot!(xᵢ, uᵢ, label="discret cos", lw=lineWidth)
end

cosAgainstDiscretCos(3)
cosAgainstDiscretCos(20)
cosAgainstDiscretCos(150)
cosAgainstDiscretCos(350)

# 16
function analogAndDiscretSolutions(n::Int)
    ω = 4 * π
    h = 1 / (n + 1)

    rangeₐ = 0:10^(-2):1

    analogSolution(x) = ((ω)^-2) * cos(ω * x) - ((ω)^-2)
    u = analogSolution.(rangeₐ)
  
    rangeᵢ = h * (1:n)
    Kₙ, Tₙ, Bₙ, Cₙ = CreateKTBC(n)
    uᵢ = h^2 * inv(Kₙ) * cos.(ω * rangeᵢ)

    plot(rangeₐ, u, label="analog", title="Analog+Discret solutions with N=$n")
    plot!(rangeᵢ, uᵢ, label="discret")
end

analogAndDiscretSolutions(4)
analogAndDiscretSolutions(8)
analogAndDiscretSolutions(16)
analogAndDiscretSolutions(32)
