"""
梯形求积公式
"""
function trapezoid(f::Function, a::Real, b::Real)
    (f(a)+f(b))*(b-a)/2
end
"""
Simpson求积公式
"""
function simpson(f::Function, a::Real, b::Real)
    c = (a+b)/2  # 中点
    (f(a)+4*f(c)+f(b))*(b-a)/3
end



f(x) = sin(x)
trapezoid(f,0,0.5)
simpson(f,0,0.5)
