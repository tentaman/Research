using DelimitedFiles
using NPZ

txt = readdlm("/home/vamuscari/Coding/Halo/test2.txt")

npzwrite("Halo200.npy", txt)

y = npzread("Halo200.npy")

println(y)