using NPZ
using Interpolations
using Polynomials, SpecialPolynomials
using Cubature
import Base.print_matrix # print_matrix()
using Plots

pyplot()
#plotly() #zoomable results


dir = homedir()
filename = "Research/camb_test.npy"
filepath = "$dir/$filename"

function main()
    camb = npzread(filepath)
    Inter = LinearInterpolation(camb[:,1], camb[:,2])
    R = LinRange( 0.0001,2.3,10000) #must be within original range  

    
    p0 = Legendre([1])
    p2 = Legendre([0, 0, 1])
    p4 = Legendre([0, 0, 0, 0, 1])
    #println("$p0 \n $p2 \n $p4")
    #println(convert(Polynomial, p2)) #will show function

    b = 0
    f = 1
    P(k, μ, pn) = Inter(k) * (b + μ^2 * f)^2 * pn(μ)

    #hcubature(x -> P(x[1], x[2], p0), [0.1,-1],[0.2,1], abstol=1e-8)

    #from 0 to 0.2 in steps of 0.01
    
    k = 0:0.01:0.2
    
    p0_bin = zeros(length(k) - 1)
    println("\n","p0 bin")
    for i in 2 :length(k) - 1
        pval = hcubature(x -> P(x[1], x[2], p0), [k[i],-1],[k[i + 1],1], abstol=1e-8)
        p0_bin[i] = pval[1] 
        println("k: ", k[i] ," to ", k[i + 1] , "  p: ", pval)
    end

    p2_bin = zeros(length(k) - 1)
    println("\n","p2 bin")
    
    for i in 2 :length(k) - 1
        pval = hcubature(x -> P(x[1], x[2], p2), [k[i],-1],[k[i + 1],1], abstol=1e-8)
        p2_bin[i] = pval[1] 
        println("k: ", k[i] ," to ", k[i + 1] , "  p2: ", pval)
    end
  
    #hcubature(x -> P(x[1], x[2], p0), [0.1,-1],[0.2,1], abstol=1e-8)

    #Plotting----------------------------------------------
    #plot(camb[:,1], camb[:,2]) #Interpolation overlap
    #display(plot!(R, Inter(R)))


    #logR = log10.(LinRange( 0.0001,2.3,10000))
    #display(plot(logR, log10.(Inter(R))))
    x = 1:20
    display(plot(x, p0_bin))
    display(plot(x, p2_bin))

    #-----------------------------------------------------



end 

main()



