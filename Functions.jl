using LinearAlgebra

# 计算三角单元的面积
function trianglearea(coordinates::Matrix{Float64})
    m, n = size(coordinates)
    @assert m == 3 && n == 2
    s = [ones(Float64, 3) coordinates]
    abs(det(s))
end
