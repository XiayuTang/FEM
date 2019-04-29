module LinearSystem  # 模板
export gauss_seidel, cg, jacobi, pcg

using LinearAlgebra

# V = Union{Vector{Int8},Vector{Int16},Vector{Int32},Vector{Int64},Vector{Int128},Vector{Float16},Vector{Float32},Vector{Float64}}
# M = Union{Matrix{Int8},Matrix{Int16},Matrix{Int32},Matrix{Int64},Matrix{Int128},Matrix{Float16},Matrix{Float32},Matrix{Float64}}
M = AbstractMatrix{<:Real}
V = AbstractVector{<:Real}

"""
    isrowzero(A)

判断矩阵A的每一行是否全为零
"""
function isrowzero(A::M)
    n = size(A,1)
    o = zeros(eltype(A),n)
    for i in 1:n
        if A[i,:] == o
            return true
        end
    end
    false
end


"""
    iscolzero(A)

判断矩阵A的每一列是否全为零
"""
function iscolzero(A::M)
    n = size(A,2)
    o = zeros(eltype(A),n)
    for i ∈ 1:n
        if A[:,i] == o
            return true
        end
    end
    false
end


"""
    diagchange!(A)

若矩阵A的对角元的某几个为0，
则通过初等行变换将0转换为非零元素
"""
function diagchange!(A::M)
    m,n = size(A)
    @assert m == n  # 要求 A为方阵
    @assert !isrowzero(A)
    @assert !iscolzero(A)
    for i ∈ 1:n
        if isapprox(A[i,i],0,atol=1e-6)
            # 若第 i个对角元素接近于 0，则寻找该对角元所处的那一列
            # 中的绝对值最大元所处的行数
            k = argmax(abs.(A[:,i]))
            if abs(A[k,i]) < 1e-6
                error("矩阵A是病态的！")
            end
            A[[i,k],:] = A[[k,i],:]  # 交换两行
        end
    end
    nothing
end


"""
    jacobi(A,b,x⁰,ϵ,maxiter) -> (x,n)

Jacobi迭代法求解线性方程组 ``Ax = b``

# Aguments
- `A`  系数矩阵
- `b`  常数项
- `x⁰`  迭代初值
- `ϵ`  允许误差
- `maxiter`  最大迭代次数
- `x`  方程组的近似解
- `n`  迭代次数
"""
function jacobi(A::M,b::V,x⁰::V,ϵ::Float64=1e-5,maxiter::Int=100) 
    diagchange!(A)
    n = 0
    D = Diagonal(diag(A))
    Dinv = inv(D)
    while true
        r⁰ = b - A*x⁰
        x = x⁰ + Dinv*r⁰
        n += 1
        if norm(x-x⁰) < ϵ
            return x, n
        end
        if n > maxiter
            println("超过最大迭代次数！")
            return nothing
        end
        x⁰ = x
    end
end


"""
    gauss_seidel(A,b,x⁰,ϵ,maxiter) -> (x,n)

Gauss-Seidel迭代法求解线性方程组 ``Ax = b``

# Aguments
- `A`  系数矩阵
- `b`  常数项
- `x⁰`  迭代初值
- `ϵ`  允许误差
- `maxiter`  最大迭代次数
- `x⁰`  方程组的近似解
- `n`  迭代次数
"""

function gauss_seidel(A::M, b::V, x⁰::V, ϵ::Float64=1e-5,maxiter::Int=100)
    diagchange!(A)
    n = 0
    L = tril(A)
    Linv = inv(L)
    while true
        r⁰ = b - A*x⁰
        x = x⁰ + Linv*r⁰
        n += 1
        if norm(x-x⁰) < ϵ
            return x, n
        end
        if n > maxiter
            println("超过最大迭代次数！")
            return nothing
        end
        x⁰ = x
    end
end
"""
    cg(A, b, x₀, e)

共轭梯度法求解 ``Ax=b``

# Aguments
- `A` 系数矩阵
- `b` 常数项
- `x₀` 迭代初值
- `e` 允许最大误差
"""
function cg(A::M, b::V, x₀::V, e::Float64=1e-5)
    @assert isposdef(A)  # 判断 A 是否为对称正定矩阵
    n = size(A, 1)
    @assert n == length(b) == length(x₀)  # 维度匹配
    d₀ = b - A*x₀
    r₀ = b - A*x₀
    for _ in 1:n
        if norm(r₀) < e
            break
        end
        α = (r₀' ⋅ r₀) / (d₀' * A * d₀)
        x = x₀ + α * d₀
        r = r₀ - α * A * d₀
        β = (r' ⋅ r) / (r₀' ⋅ r₀)
        d = r + β * d₀
        r₀ = r
        x₀ = x
        d₀ = d
    end
    x₀
end
"""
    pcg(A, b, x₀, e)

预处理共轭梯度法求解 ``Ax=b``

# Aguments
- `A` 系数矩阵
- `b` 常数项
- `x₀` 迭代初值
- `e` 允许最大误差
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


end  # end Module