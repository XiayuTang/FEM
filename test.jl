# test codes and functions here
# how to use debugger on Juno
# step1: make sure there is a function called main() which runs your
# entire program.
# step2: using Revise, Juno
# step3: includet("Your File Name")
# step3: Juno.@enter main() or Juno.@run main()
# push!(LOAD_PATH,".");
using LinearSystem
# function main()
#     a = 10
#     b = 20
#     A = [1 3 4 5]
#     B = rand(4,4)
#     α = a*A
#     β = b*B
#     α * β
#     nothing
# end
a = [1.0 0.0; 0.0 1.0];
b = [1.0;2.0];
pcg(a,b,zeros(2))