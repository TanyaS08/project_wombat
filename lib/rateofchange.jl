"""
    _rate_gradient(∂𝑋, ∂𝑌)

Returns the rate of change in units of the values stored in the
grid, and the angle of the change in degrees. When both ∂X and ∂Y are
equal to 0, the angle is assumed to be 0.
"""
function _rate_gradient(∂𝑋::T, ∂𝑌::T) where {T <: Number}
    if ∂𝑋 == ∂𝑌 == 0.0
        return (0.0, 0.0)
    end
    m = sqrt(∂𝑋^2 + ∂𝑌^2)
    Δ = ∂𝑋 < 0.0 ? 180.0 : 0.0
    θ = rad2deg(atan(∂𝑋 / ∂𝑌)) + Δ
    θ = isnan(θ) ? 0.0 : θ
    return (m, θ)
end

"""
    _rateofchange(x::Vector{T}, y::Vector{T}, z::Vector)

Rate of change for a series of three points, defined as a series of `x` and `y`
coordinates and a value `z`. Returns a rate of change (in units of `z`) and a
gradient in degrees.
"""
function _rateofchange(x::Vector{T}, y::Vector{T}, z::Vector) where {T<:Number}

    # Check that all three vectors have the same length
    length(x) == length(y) || throw(DimensionMismatch("x and y must have the same length"))
    length(x) == length(z) || throw(DimensionMismatch("x and z must have the same length"))

    # Get the matrix of coefficients
    C = cat(y, x, fill(one(T), length(x)); dims=(2, 2))
    coeff = Base.inv(C) * z

    𝑋 = sum(C[:, 1]) / 3.0
    𝑌 = sum(C[:, 2]) / 3.0

    ∂𝑋 = coeff[2] * 𝑌 + coeff[3]
    ∂𝑌 = coeff[1] * 𝑋 + coeff[3]

    # Rate of change and direction
    return _rate_gradient(∂𝑋, ∂𝑌)
end

"""
    _rateofchange(A::Matrix{T}) where {T <: Number}

Returns the rate of change and the gradient for a 2x2 grid of numbers.
"""
function _rateofchange(A::Matrix{T}) where {T<:Number}
    size(A) == (2, 2) || throw(DimensionMismatch("the matrix A must have size (2,2)"))

    # We can get the values directly from the matrix
    𝑍₄, 𝑍₁, 𝑍₃, 𝑍₂ = A

    ∂𝑋 = 𝑍₂ - 𝑍₁ + 0.5(𝑍₁ - 𝑍₂ + 𝑍₃ - 𝑍₄)
    ∂𝑌 = 𝑍₄ - 𝑍₁ + 0.5(𝑍₁ - 𝑍₂ + 𝑍₃ - 𝑍₄)

    # Rate of change and direction
    return _rate_gradient(∂𝑋, ∂𝑌)
end