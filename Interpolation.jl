using NPZ
using Interpolations
using Polynomials, SpecialPolynomials
import Base.print_matrix # print_matrix()
using Plots
pyplot()
#plotly() #zoomable results



function main()
    camb = npzread("/home/van/Research/camb_test.npy")
    Inter = LinearInterpolation(camb[:,1], camb[:,2])
    R = LinRange( 0.0001,2.3,10000) #must be within original range  

    #P(k) = Inter(k) * (b * μ^2 * f)^2
    #P(k) = Inter(k) * (b * μ^2 * f)^2

    p0 = Legendre([1])
    p2 = Legendre([0, 0, 1])
    p4 = Legendre([0, 0, 0, 0, 1])
    println("$p0 \n $p2 \n $p4")
    #println(convert(Polynomial, p2)) #will show function
    




    #Plotting----------------------------------------------
    plot(camb[:,1], camb[:,2]) #Interpolation overlap
    display(plot!(R, Inter(R)))


    logR = log10.(LinRange( 0.0001,2.3,10000))
    display(plot(logR, log10.(Inter(R))))
    #-----------------------------------------------------





end 

main()



