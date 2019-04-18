using LinearAlgebra
function trianglearea(coordinates::Matrix{Float64})
    m, n = size(coordinates)
    @assert m == 3 && n == 2
    s = [ones(Float64, 3) coordinates]
    abs(det(s))
end
