module SpecialVectors
using Test

function constants(n=5)
    return ones(n)
end

@test constants() == ones(5)
@test constants(6) == ones(6)

function linear(n=5)
    return [1:n...]
end

@test linear() == [1:5...]
@test linear(6) == [1:6...]

function squares(n=5)
    return [1:n...].^2
end

@test squares() == [1:5...].^2
@test squares(6) == [1:6...].^2

function checkDim(k::Int, n::Int)
    if (k > n)
        throw(DimensionMismatch("k must be less than n!"))
    end
    return true
end

@test_throws DimensionMismatch checkDim(6, 1)
@test checkDim(1, 5)

function deltaAtK(k=1, n=5)
    checkDim(k, n)
    return [zeros(k - 1)..., 1, zeros(n - k)...]
end

@test_throws DimensionMismatch deltaAtK(6)
@test deltaAtK() == [1, 0, 0, 0, 0]
@test deltaAtK(3, 6) == [0, 0, 1, 0, 0, 0]

function stepAtK(k=1, n=5)
    checkDim(k, n)
    return [zeros(k - 1)..., 1, ones(n - k)...]
end

@test_throws DimensionMismatch stepAtK(6)
@test stepAtK() == [1, 1, 1, 1, 1]
@test stepAtK(3) == [0, 0, 1, 1, 1]

function rampAtK(k=1, n=5)
    checkDim(k, n)
    return [zeros(k - 1)..., 0:(n - k)...]
end

@test_throws DimensionMismatch rampAtK(6)
@test rampAtK() == [0, 1, 2, 3, 4]
@test rampAtK(2) == [0, 0, 1, 2, 3]
@test rampAtK(3) == [0, 0, 0, 1, 2]

function sines(t, n=5)
    return sin.([1:n...] .* t)
end

@test round.(sines(π / 2), digits=2) == [1, 0, -1, 0, 1]
@test round.(sines(π / 4), digits=2) == [.71, 1.0, .71, 0, -.71]

function cosines(t, n=5)
    return cos.([1:n...] .* t)
end

@test round.(cosines(π / 2), digits=2) == [0, -1, 0, 1, 0]
@test round.(cosines(π / 4), digits=2) == [.71, 0.0, -.71, -1.0, -.71]

function exponentials(t, n=5)
    return exp.([1:n...] .* (im * t))
end

@test exponentials(0) == repeat([1.0 + 0.0im], 5)
@test exponentials(1) == [
  0.5403023058681398 + 0.8414709848078965im
 -0.4161468365471424 + 0.9092974268256817im
 -0.9899924966004454 + 0.1411200080598672im
 -0.6536436208636119 - 0.7568024953079282im
 0.28366218546322625 - 0.9589242746631385im
]

export constants, linear, squares, deltaAtK, stepAtK, rampAtK, sines, cosines, exponentials

end