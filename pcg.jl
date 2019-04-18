using LinearAlgebra
function pcg(A::Matrix{T}, b::Vector{T}, x₀::Vector{T},e::Float64=1e-5) where T <: AbstractFloat
    @assert isposdef(A)
    n = size(A, 1)
    @assert n == length(b) == length(x₀)
    M = Diagonal(diag(A))
    r₀ = b - A*x₀
    Minv = inv(M)
    d₀ = Minv*r₀
    z₀ = d₀
    for _ in 1:n
        if norm(r₀) < e
            break
        end
        α = (r₀' ⋅ z₀) / (d₀' * A * d₀)
        x = x₀ + α * d₀
        r = r₀ - α*A*d₀
        z = Minv * r
        β = (r'⋅z)/(r₀'⋅z₀)
        d = z + β*d₀
        x₀ = x
        r₀ = r
        z₀ = z
        d₀ = d
    end
    x₀
end  # function pcg
