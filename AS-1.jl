##
x, y, z = symbols("x y z")
(x, y, z)
solveset(x^2-1-y, x)
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

##
let 
    π, ω1A, ω1B, ω2A, ω2B = 0.5, 2, 5, 5, 2
    solveset(π * (ω1A + p*ω1B) + ((ω2A + p*ω2B) * π ) - ω1A - ω2A , p)
    end 
##

1/3
