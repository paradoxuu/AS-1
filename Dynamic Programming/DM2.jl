using Plots
using BenchmarkTools

function vfsolvex(vnew, kgrid, tolerance, imax, σ)
    β = 0.9

    v = vnew .+ 2*tolerance
    cartesianindex = Array{CartesianIndex{2}, length(v)}
    i = 1
    
    c = kgrid .- kgrid'
    c[c .< 0] .= 0
    u = (c .^ (1-σ) .- 1) ./ (1-σ)
    u[kgrid .- kgrid' .< 0] .= -Inf

    while maximum(abs.(v - vnew)) > tolerance && i <= imax
        v = vnew;
        
        (vnew, cartesianindex) = findmax(u .+ β*v', dims = 2);
        i += 1;
    end
    kprimeindex = getindex.(cartesianindex, 2)
    return (v = vnew, kprime = kgrid[kprimeindex], kprimeindex = kprimeindex)
end


kupper = 5
klower = 0.01
n = 1000
kgrid = collect(range(klower, stop = kupper, length = n))
(v, kprime, kprimeindex) = vfsolvex(zeros(n), kgrid, 0.001, 1000, 1.5);

scatter(kgrid, v, label = "v")
scatter(kgrid, kprime, label = "k'")

function policy(x)
    e=getindex(findall(kgrid .== x),1)
    return kprime[e]  
end 

T = 100
kpath = zeros(T)
k = kgrid[1000]
for i in 1:T
    k = policy(k)
    kpath[i] = k
end 
time = collect(range(1, T, length = T))
scatter(time, kpath, label = "kpath")
