α = 0.3
δ = 0.1
β = 0.9
k_l = 2
k_u = 0.001
n = 10
kgrid = collect(range(k_l, stop = k_u, length =n))

v = zeros(n)
tolerance = 0.001
imax = 1000
iter = 1
vnew = zeros(n)
while maximum(abs.(v-vnew) > tolerance) && iter < imax
v=vnew
c = zeros(n,n) 
for i in 1:n
    for j in 1:n
    c[i,j] = kgrid[i]^α + (1-δ)*kgrid[i] - kgrid[j]
        if c[i,j]<0
            c[i,j]=0
        end 
    end
end 
vnew = maximum(log.(c) .+ β*v' , dim = 2)
iter += 1
end 
