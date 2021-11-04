using Plots
using BenchmarkTools

function vfsolve(vnew, kgrid, tolerance, imax)
    α = 0.37
    δ = 0.1
    β = 0.95
    σ = 1.5

    v = vnew .+ 2*tolerance
    cartesianindex = Array{CartesianIndex{2}, length(v)}
    i = 1
    
    c = kgrid.^α + (1-δ)*kgrid .- kgrid'
    c[c .< 0] .= 0
    u = (c .^ (1-σ) .- 1) ./ (1-σ)
    u[kgrid.^α + (1-δ)*kgrid .- kgrid' .< 0] .= -Inf

    while maximum(abs.(v - vnew)) > tolerance && i <= imax
        v = vnew;
        
        (vnew, cartesianindex) = findmax(u .+ β*v', dims = 2);
        i += 1;
    end
    kprimeindex = getindex.(cartesianindex, 2)
    return (v = vnew, kprime = kgrid[kprimeindex], kprimeindex = kprimeindex)
end


kupper = 5
klower = 0.05
n = 1000
kgrid = collect(range(klower, stop = kupper, length = n))
(v, kprime, kprimeindex) = vfsolve(zeros(n), kgrid, 0.001, 1000);

plot(kgrid, kprime, label = "k'")


function policy(x)
    e=getindex(findall(kgrid .== x),1)
    return kprime[e]  
end 

T = 100
kpath = zeros(T)
k = kgrid[92]
for i in 1:T
    k = policy(k)
    kpath[i] = k
end 
time = collect(range(1, T, length = T))
plot(time, kpath, label = "kpath")



# Check k bounds, stepsize, tolerance, imax
findall(kprime .== maximum(kgrid))
any(kprime .== minimum(kgrid))




function vfsolve1(vnew, kgrid, imax)
    α = 1.5
    δ = 0.1
    β = 0.95
    σ = 1.5

    cartesianindex = Array{CartesianIndex{2}, length(v)}
    i = 1
    
    c = kgrid.^α + (1-δ)*kgrid .- kgrid'
    c[c .< 0] .= 0
    u = (c .^ (1-σ) .- 1) ./ (1-σ)
    u[kgrid.^α + (1-δ)*kgrid .- kgrid' .< 0] .= -Inf

    while i <= imax

        (vnew, cartesianindex) = findmax(u .+ β*vnew', dims = 2);
        i += 1;

    end
    kprimeindex = getindex.(cartesianindex, 2)
    return (v = vnew, kprime = kgrid[kprimeindex], kprimeindex = kprimeindex)
end

(v1, kprime1, kprimeindex1) = vfsolve1(zeros(n), kgrid, 1);
(v5, kprime5, kprimeindex5) = vfsolve1(zeros(n), kgrid, 5);
(v10, kprime10, kprimeindex10) = vfsolve1(zeros(n), kgrid, 10);
scatter(kgrid, v1, label = "v1")
scatter(kgrid, v5, label = "v5")
scatter(kgrid, v10, label = "v10")



function vfsolve2(vnew, kgrid, tolerance, imax)
    α = 1.5
    δ = 0.1
    β = 0.95
    σ = 0.5

    v = vnew .+ 2*tolerance
    cartesianindex = Array{CartesianIndex{2}, length(v)}
    i = 1
    
    c = kgrid.^α + (1-δ)*kgrid .- kgrid'
    c[c .< 0] .= 0
    u = (c .^ (1-σ) .- 1) ./ (1-σ)
    u[kgrid.^α + (1-δ)*kgrid .- kgrid' .< 0] .= -Inf

    while maximum(abs.(v - vnew)) > tolerance && i <= imax
        v = vnew;
        
        (vnew, cartesianindex) = findmax(u .+ β*v', dims = 2);
        i += 1;
    end
    kprimeindex = getindex.(cartesianindex, 2)
    return (v = vnew, kprime = kgrid[kprimeindex], kprimeindex = kprimeindex)
end

(v, kprime, kprimeindex) = vfsolve2(zeros(n), kgrid, 0.001, 1000);
scatter(kgrid, v, label = "v")
scatter(kgrid, kprime, label = "k'")