from __future__ import division
from math import sqrt
from math import exp
from math import pow as p
import numpy as np
import matplotlib.pyplot as plt
import time

V = np.zeros(4)
V_new = np.zeros(4)

K = 0.971635

tp = 2*np.pi

steps = 1e3

m = 1
n = -1

linIC_density = 60
lICd = np.array([20, 40, 60, 80, 100, 200])
initials = int(p(linIC_density, 2))

initConds = np.zeros([initials, 2])

for j in range(0, initials):
    initConds[j, 0] = (j/linIC_density)%1.0
    initConds[j, 1] = np.floor(j/linIC_density)/linIC_density

def standMap(v):
    V_new[1] = v[1] + K*np.sin(tp*v[0])/tp
    V_new[0] = v[0] + V_new[1]

    V_new[3] = v[2]*K*np.cos(tp*v[0])/tp + v[3]
    V_new[2] = V_new[3] + v[2]
    
    return V_new

renorm = 10
k = 0

Y = np.zeros(int(steps/renorm))
Y_ave = np.zeros(int(steps/renorm))
Y[0] = 1
Y_ave[0] = 1

t0 = time.time()
i = 1

chaotic = 0
times = np.linspace(0, steps, int(steps/renorm))
A = np.vstack([times, np.ones(len(times))]).T
a = 0
b = 0

numPoints = int(np.size(lICd))
prevalence = np.zeros(numPoints)





for q in range(0, numPoints):
    linIC_density = lICd[q]
    chaotic = 0
    initials = int(p(linIC_density, 2))

    initConds = np.zeros([initials, 2])

    for j in range(0, initials):
        initConds[j, 0] = (tp*j/linIC_density)%tp
        initConds[j, 1] = tp*np.floor(j/linIC_density)/linIC_density
        
    for l in range(0, initials):
        dev = np.array([np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None)])
        mag = sqrt(p(dev[0], 2) + p(dev[1], 2))
        dev = dev/mag
        mag = sqrt(p(dev[0], 2) + p(dev[1], 2))
        
        V = np.concatenate((initConds[l, :], dev))
        
        i = 1
        while i < int(steps/renorm):
            k = 0
            while k < renorm:
                V = standMap(V)
            
                V[0] = V[0]%1.0
                V[1] = V[1]%1.0
        
                k += 1
        
            dev = np.array([V[2], V[3]])
            mag = sqrt(p(dev[0], 2) + p(dev[1], 2))
        
            Y[i] = Y[i - 1]/(m + 1)/p(i, n)
            Y[i] = Y[i] + np.log(mag)*p(i + 1, m)
            Y[i] = (m + 1)*p(i + 1, n)*Y[i]
            dev = dev/mag
        
            Y_ave[i] = p(i, m + n + 1)*Y_ave[i - 1]
            Y_ave[i] = Y_ave[i] + Y[i]
            Y_ave[i] = Y_ave[i]/p(i + 1, m + n + 1)
        
            V[2] = dev[0]
            V[3] = dev[1]
            
            i += 1
            
        a, b = np.linalg.lstsq(A, 2*Y_ave, rcond = None)[0]
        a = np.log10(np.abs(a))
        
        if a > -2.2:
            chaotic += 1
            
    prevalence[q] = chaotic/initials
    print 'Prop: ', prevalence[q]
    t1 = time.time()
    print 'Time: ', t1 - t0
    print





t1 = time.time()
print 'Total time Elapsed: ', t1 - t0

print(prevalence)
plt.plot(lICd, prevalence)
plt.xlabel('Linear Density of ICs')
plt.ylabel('Estimated Proportion of Chaotic Orbits')

plt.show()
