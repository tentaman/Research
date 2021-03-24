#For reading Numpy array, for error function, for random number.
using NPZ
using SpecialFunctions
using Random: rand
using Distributions
#corfunc

#Calling Halos by length of array because it will have to itterate through the coordinates as well. 
#Thought about breaking up main into more functions but I'm doubtful it will provide any additional performance


function main()
    Halo   = npzread("Halo200.npy")
    len    = length(Halo[:,1])
    Total  = 0
    Mmin   = 10e12
    Mprime = 10e13
    σ      = 0.15 
    α      = .93
    #Galaxy Array
    #Mass of galaxy, position x, position y, position z, Central or satellite
    Galaxy = ["Mass of galaxy" "position x" "position y" "position z" "Central or satellite"]


    for i in range(1, length = len)
        CentralP = CentralProbability(Halo[i,1], Mmin, σ)
        if CentralP >= rand(Float64)
            CentGal    = [Halo[i,1] Halo[i,4] Halo[i,5] Halo[i,6] "C"]
            Galaxy     = vcat(Galaxy, CentGal)
            SatelliteP = SatelliteProbability(Halo[i,1], Mmin, Mprime, α)
            if SatelliteP >= 1
                for s in range(1, length = SatelliteP)
                    SXYZ   = PointGen(Halo[i,2], Halo[i,4], Halo[i,5], Halo[i,6])
                    SatGal = [10e10 SXYZ "S"]
                    Galaxy = vcat(Galaxy, SatGal)
                end
            end
        end

    end
    
    show(stdout, "text/plain", Galaxy)

end

#will generate random number between (-1,1)
#a sort of workaround for now
function RandNum()
    num = rand()
    if rand() < 0.5
        num = -num
    end
    return num
end

function CentralProbability(M, Mmin, σ)
    p = (1/2) * (1 + erf((log10(M)-log10(Mmin))/σ))
    return p
end

function SatelliteProbability(M, Mmin, Mprime, alpha)
    #Since Σ P(x) = 1, will add poisson until it surpasses rand()
    #poisson has to start at P(X = 0) to be inclusive with while loop
    MeanAve = ((M - Mmin) / Mprime)^alpha
    SatNum  = 0
    poisson = pdf(Poisson(MeanAve),0)
    check   = rand()
    
    while poisson < check
        poisson += pdf(Poisson(MeanAve),SatNum) 
        SatNum  += 1   
    end

    return SatNum
end

function PointGen(r, x, y, z)
    #unit vector generation
    ux = 0
    uy = 0
    uz = 0
    while true 
        ux = RandNum()
        uy = RandNum()
        uz = RandNum()
        if ( sqrt( (ux^2) + (uy^2) + (uz^2) ) < 1 )
            break
        end
    end

    #translation of coordinates
    x = r*ux + x
    y = r*uy + y
    z = r*uz + z

    return [x y z]

end

main()
println("\n")