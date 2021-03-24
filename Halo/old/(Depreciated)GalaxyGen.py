import numpy as np
import random
import math
import time

t0 = time.time()

txt = np.loadtxt('/home/vamuscari/Coding/Halo/test2.txt')
t1 = time.time()

print("Load Time: " + str((t1 - t0)/60))

Mass = [int(m) for m in txt[:,0]]
MassTF = []
Total = 0


Mmin = 10 * (10 ** 12)
σ = 0.15 

for M in Mass:
    Num = (1/2)* (1 + math.erf((np.log10(M)-np.log10(Mmin))/σ))
    if Num >= random.random():
        MassTF.append([M,1])
        Total += 1
    else:
        MassTF.append([M,0])

print(Total)