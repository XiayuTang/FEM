# test codes and functions here
# how to use debugger on Juno
# step1: make sure there is a function called main() which runs your
# entire program.
# step2: using Revise, Juno
# step3: includet("Your File Name")
# step3: Juno.@enter main() or Juno.@run main()
function main()
    a = 10
    b = 20
    A = [1 3 4 5]
    B = rand(4,4)
    α = a*A
    β = b*B
    α * β
end
