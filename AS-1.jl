##
x, y, z = symbols("x y z")
(x, y, z)
k = solveset(x^4-x^3+5*x^2+1, x)
##
##
a, b, c = symbols("a b c")
nonlinsolve([a^2 + a, a - b - c], a, b)
##
π, p, γ, ω1A, ω1B, ω2A, ω2B = symbols("π p γ ω1A ω1B ω2A ω2B")
solveset(π * (ω1A + p*ω1B) + ((ω2A + p*ω2B) / (1 + p * ((1 - π)/(p * π))^(1/γ))) - ω1A - ω2A , p)
##
let 
π, γ, ω1A, ω1B, ω2A, ω2B = 0.5, 1, 2, 5, 5, 2
solveset(π * (ω1A + p * ω1B) + (ω2A + p * ω2B) / (1 + p * ((1 - π)/(p * π))^(1/γ)) - ω1A - ω2A , p)
end
##
us
##
let 
    π, ω1A, ω1B, ω2A, ω2B = 0.5, 2, 5, 5, 2
    solveset(π * (ω1A + p*ω1B) + ((ω2A + p*ω2B) * π ) - ω1A - ω2A , p)
    end 
##

##
let 
    π, γ, ω1A, ω1B, ω2A, ω2B = 0.5, 2, 2, 5, 5, 2
    solveset(π * (ω1A + p * ω1B) + (ω2A + p * ω2B) / (1 + p * ((1 - π)/(p * π))^(1/γ)) - ω1A - ω2A , p)
    end
##

##
let 
    π, γ, ω1A, ω1B, ω2A, ω2B = 0.5, 0.5, 2, 5, 5, 2
    solveset(π * (ω1A + p * ω1B) * (1 + p * ((1 - π)/(p * π))^(1/γ)) + (ω2A + p * ω2B) - (ω1A + ω2A) * (1 + p * ((1 - π)/(p * π))^(γ)) , p)
    end
##

##
using NLsolve

function f!(F, x)
    F[1] = (x[1]+3)*(x[2]^3-7)+18
    F[2] = sin(x[2]*exp(x[1])-1)
end

function j!(J, x)
    J[1, 1] = x[2]^3-7
    J[1, 2] = 3*x[2]^2*(x[1]+3)
    u = exp(x[1])*cos(x[2]*exp(x[1])-1)
    J[2, 1] = x[2]*u
    J[2, 2] = u
end

nlsolve(f!, j!, [ 0.1; 1.2])
##
##
using LineSearches, Optim

ϕ(x) = (x - π)^4
dϕ(x) = 4*(x-π)^3
ϕdϕ(x) = ϕ(x),dϕ(x)

α0 = 9.0
ϕ0 = ϕ(0.0)
dϕ0 = dϕ(0.0)

for ls in (Static,BackTracking,MoreThuente,StrongWolfe)
    println(ls, ": ", (ls())(ϕ, dϕ, ϕdϕ, α0, ϕ0,dϕ0))
end

##
##
using Roots
let 
    π, γ, ω1A, ω1B, ω2A, ω2B = 0.5, 1, 2, 5, 5, 2
    f(p) = π * (ω1A + p*ω1B) + ((ω2A + p*ω2B) / (1 + p * ((1 - π)/(p * π))^(1/γ))) - ω1A - ω2A
    find_zero(f, 1.5)
end 
##

function GE(π, γ, ω1A, ω1B, ω2A, ω2B)
    f(p) = π * (ω1A + p*ω1B) + ((ω2A + p*ω2B) / (1 + p * ((1 - π)/(p * π))^(1/γ))) - ω1A - ω2A
    p = round(find_zero(f, 1.5), digits= 3)
    C1A = π * (ω1A + p*ω1B)
    C1B = ((1 - π) / p) * (ω1A + p*ω1B)
    C2A = (ω2A + p*ω2B) / (1 + p * ((1 - π)/(p * π))^(1/γ))
    C2B = ((ω2A + p*ω2B) / (1 + p * ((1 - π)/(p * π))^(1/γ))) * ((1 - π)/(p * π))^(1/γ)
    @show p ,C1A ,C1B, C2A, C2B   
end 
##
#4(e) No aggregate uncertanity 
GE(0.5, 2, 2, 5, 5, 2)
#4(f) There is no aggregate uncertainty but for person A the good world is more likely to happen 
GE(0.25, 2, 2, 5, 5, 2)
#4(f) There is no aggregate uncertainty but for person B the good world is more likely to happen
GE(0.75, 2, 2, 5, 5, 2)
#4(g)(e) γ=1.5
GE(0.5, 1.5, 2, 5, 5, 2)
#4(g)(e) γ=3
GE(0.5, 3, 2, 5, 5, 2)
#4(g)(f) γ=1.5 π=0.25
GE(0.25, 1.5, 2, 5, 5, 2)
#4(g)(f) γ=1.5 π=0.75
GE(0.75, 1.5, 2, 5, 5, 2)
#4(g)(f) γ=3 π=0.25
GE(0.25, 3, 2, 5, 5, 2)
#4(g)(f) γ=3 π=0.75
GE(0.75, 3, 2, 5, 5, 2)
#4(g)(e) γ=2 π=0.5
GE(0.5, 2, 2, 3, 5, 2)
#4(g)(e) γ=1.5 π=0.5
GE(0.5, 1.5, 2, 3, 5, 2)
#4(g)(e) γ=3 π=0.5
GE(0.5, 3, 2, 3, 5, 2)


