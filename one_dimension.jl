# push!(LOAD_PATH,".")
# using Integral 
# using LinearSystem
using LinearAlgebra
M = AbstractMatrix{<:Real}
V = AbstractVector{<:Real}
"""
预处理共轭梯度法
"""
function pcg(A::M, b::V, x₀::V,e::Float64=1e-5)
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
end

"""
Simpson求积公式
"""
function simpson(f::Function, a::Real, b::Real)
    c = (a+b)/2  # 中点
    (f(a)+4*f(c)+f(b))*(b-a)/3
end

"""
网格（区间）划分
"""
function mesh(a::Real, b::Real, h::Float64=0.01)
    return collect(a:h:b)
end

"""
求出刚度矩阵
"""
function stiff(p::Function, q::Function, xmesh, h)
    n = length(xmesh) - 1
    a₁ = zeros(Float64, n-1)
    a₂ = zeros(Float64, n)
    for i in 1:n-1
        f(x) = -p(xmesh[i]+h*x)/h + h*q(xmesh[i]+h*x)*x*(1-x) 
        a₁[i] = simpson(f,0,1)
    end
    for k in 1:n
        g(x) = p(xmesh[k]+h*x)/h + h*q(xmesh[k]+h*x)*x^2 + p(xmesh[k+1]+h*x)/h + h*q(xmesh[k+1]+h*x)*(1-x)^2
        a₂[k] = simpson(g,0,1)
    end
    A = diagm(0=>a₂, 1=>a₁, -1=>a₁)
    A
end
