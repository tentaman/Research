#For reading Textfile, for error function, for random number.
using NPZ
using SpecialFunctions
using Random

#use package NPZ for reading writing array
#julia uses jld package 

@time begin


txt = npzread("Halo200.npy")

Mass = [m for m in txt[:,1]]
MassTF = [0 0]
Total = 0


Mmin = 10 * (10 ^ 12)
σ = 0.15 



for M in Mass
    Num = (1/2) * (1 + SpecialFunctions.erf((log10(M)-log10(Mmin))/σ))
    if Num >= rand(Float64)
        global MassTF = vcat(MassTF, [M 1])
        global Total += 1
            
    else
        global MassTF = vcat(MassTF, [M 0])
    end
end
println(Total)

end