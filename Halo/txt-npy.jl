#-------------------------------------------------------------------------------
#txt to npy writer-
#or npy to txt writer.
#Written using ARGS[] so this script is meant to be used from the command line.
#Pass file into script from command line.
#File must be in same folder as script or must be called from directory.
#Examples:
#julia txt-npy.jl data.txt
#julia txt-npy.jl data.npy
#julia /home/[username]/txt-npy.jl /home/[username]/data.txt
#-------------------------------------------------------------------------------

using DelimitedFiles
using NPZ


function main()
    bool , type = Errorcheck() 
    if bool  == true
        println("Writing...")
        @time Write(type)
    end
end

#Read file then convert file
function Write(type)
    l = length(ARGS[1])
    name = SubString(ARGS[1], 1:(l-4))
    if  type == "txt"
        txt = readdlm(ARGS[1])
        name  = join([name , ".npy"])
        println(name)
        npzwrite(name, txt)

    elseif type == "npy"
        npy = npzread(ARGS[1])
        name  = join([name , ".txt"])
        open(name, "w") do io
            writedlm(io, npy)
        end
    end
end

#Check for correct file types
function Errorcheck()
    l = length(ARGS[1])
    type = SubString(ARGS[1], (l-2):l)
    if type != "npy" && type != "txt"
        println("Incorrect file type")
        println("Only accepts .txt and .npy")
        return false
    else
        return true , type
    end
end

main()