module DifferencesMatrix

using Test

import LinearAlgebra.diagm

function backwardDiff(dim::Int)
    return diagm(-1 => -ones(dim - 1), 0 => ones(dim))
end

@test backwardDiff(3) == [
  1 0 0
  -1 1 0
  0 -1 1
]

function forwardDiff(dim::Int)
    return diagm(0 => -ones(dim), 1 => ones(dim - 1))
end

@test forwardDiff(3) == [
  -1 1 0
  0 -1 1
  0 0 -1
]

@test backwardDiff(3) * forwardDiff(3) == [
  -1 1 0
  1 -2 1
  0 1 -2
]

function centeredDiff(dim::Int)
    return diagm(-1 => -ones(dim - 1), 0 => zeros(dim), 1 => ones(dim - 1))
end

@test centeredDiff(3) == backwardDiff(3) - (-forwardDiff(3))

export backwardDiff, forwardDiff, centeredDiff

end
