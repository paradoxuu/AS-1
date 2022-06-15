##
using Roots
function GE(π, γ, ω1A, ω1B, ω2A, ω2B)
    f(p) = π * (ω1A + p*ω1B) + ((ω2A + p*ω2B) / (1 + p * ((1 - π)/(p * π))^(1/γ))) - ω1A - ω2A
    p = round(find_zero(f, 1.5), digits= 3)
    C1A = π * (ω1A + p*ω1B)
    C1B = ((1 - π) / p) * (ω1A + p*ω1B)
    C2A = (ω2A + p*ω2B) / (1 + p * ((1 - π)/(p * π))^(1/γ))
    C2B = ((ω2A + p*ω2B) / (1 + p * ((1 - π)/(p * π))^(1/γ))) * ((1 - π)/(p * π))^(1/γ) 
    print("p = $p \n", "C1A = $C1A \n", "C1B = $C1B \n", "C2A = $C2A \n", "C2B = $C2B \n")
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
#hahahah