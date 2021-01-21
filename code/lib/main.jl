"""
    Description
"""
function LinearAlgebra.rank(N::T) where {T <: DeterministicNetwork}
    return rank(N.A)
end
