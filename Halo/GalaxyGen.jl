using NPZ                #read write
using SpecialFunctions   #error function
using Random: rand       #Random
using Distributions      #Poisson Distribution
using PyCall             #Calling python packages
import Base.print_matrix #Printing less gross matrix


#Corrfunc
np = pyimport("numpy")                  
theory = pyimport("Corrfunc.theory")


#Calling Halos by length of array because it will have to itterate through the coordinates as well. 


function main()
    Halo   = npzread("/home/van/halo/Halo500.npy")
    len    = length(Halo[:,1])
    print("length: $len \n")
    Mmin   = 10^12.69
    M0     = 10^12.94
    Mprime = 10^13.82
    σ      = 0.15
    α      = 1.08
    Sattotal = 0
    #Galaxy Array
    #Mass of galaxy, position x, position y, position z, Central or satellite
    Galaxy = [0 0 0]

    for i in range(1, length = len)
        CentralP = CentralProbability(Halo[i,1], Mmin, σ)
        if CentralP >= rand(Float64)
            CentGal    = [Halo[i,4] Halo[i,5] Halo[i,6]]
            Galaxy     = vcat(Galaxy, CentGal)
            SatelliteP = SatelliteProbability( Halo[i,1] , M0, Mprime, α)
            if SatelliteP >= 1
                for s in 1:SatelliteP
                    SXYZ   = PointGen(Halo[i,2], Halo[i,4], Halo[i,5], Halo[i,6])
                    Sattotal += 1
                    Galaxy = vcat(Galaxy, round.(SXYZ, digits=6))
                end
            end
        end

    end
    Boxsize = BoxFinder(Galaxy)
    XI(Galaxy, Boxsize)
    WP(Galaxy, Boxsize)
    #print("Number of Galaxies: $(length(Galaxy[:,1])) \n")
    #print("Satalites: $Sattotal \n")
    #print_matrix(stdout, Galaxy)

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

function SatelliteProbability(M, M0, Mprime, alpha)
    #Since Σ P(x) = 1, will add poisson until it surpasses rand()
    #poisson has to start at P(X = 0) to be inclusive with while loop
    if M <= M0
        return 0
    end
    MeanAve = ((M - M0) / Mprime)^alpha
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
    r = r/1000
    x = r*ux + x
    y = r*uy + y
    z = r*uz + z

    return [x y z]

end


#Find maximim radius in Galaxy list for Box size in corrfunc
function BoxFinder(Galaxy)
    len2 = length( Galaxy[:,1] )
    r = zeros(len2)
    for i in 1:len2
        r[i] = sqrt(sum(Galaxy[i,1:2].^2))
    end
    return (maximum(r))
end



function XI(Galaxy, Boxsize)
    nthreads = 12
    rbins = np.logspace(np.log10(0.1), np.log10(10.0), 15)
    X = Galaxy[:,1]
    Y = Galaxy[:,2]
    Z = abs.(Galaxy[:,3])
    results = theory.xi(Boxsize, nthreads, rbins, X, Y, Z)
    return results
end 


function WP(Galaxy, Boxsize)
    nthreads = 12
    rbins = np.logspace(np.log10(0.1), np.log10(10.0), 15)
    pimax = 1
    X = Galaxy[:,1]
    Y = Galaxy[:,2]
    Z = abs.(Galaxy[:,3])
    results = theory.wp(Boxsize, pimax, nthreads, rbins, X, Y, Z)
    return results
end 


main()
println("\n")


#--------------------------------------------
#   Information of mdhalo_z038_mrrs_xyz_vxyz
#   length of array: 162854827 
#   Maximum Radius: 4324.448891506933
#   Maximum Mass: 3.174e15
#   Minimum Mass: 4.7044e10
#-------------------------------------------