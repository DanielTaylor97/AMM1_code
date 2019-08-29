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

mGrid = np.zeros(initials)
times = np.linspace(0, steps, int(steps/renorm))

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

    A = np.vstack([times, np.ones(len(times))]).T
    a, b = np.linalg.lstsq(A, 2*Y_ave, rcond = None)[0]
    mGrid[l] = a

    if l%10 == 0:
        print((l + 1)/initials)

        


t1 = time.time()
print 'Time Elapsed: ', t1 - t0




#f= open("2DSMcontour.txt","w+")
#
#for k in range(0, initials):
#    f.write(str(mGrid[k]))
#    f.write("\n")
#
#f.close()




grid = np.zeros((linIC_density, linIC_density))
chaotic = 0
total = p(linIC_density, 2)

for x in range(0, linIC_density):
    for y in range(0, linIC_density):
        grid[x, y] = np.log10(np.abs(mGrid[x*linIC_density + y]))

        if grid[x, y] >-2.2:
            chaotic += 1

print 'Proportion of orbits that are chaotic: ', chaotic/total

x1 = np.linspace(0, 1.0, linIC_density)
x2 = np.linspace(0, 1.0, linIC_density)

N = 100
plt.contourf(x1, x2, grid, N, cmap = "viridis")
plt.colorbar(orientation='horizontal')
plt.xlabel('$x_1$')
plt.ylabel('$x_2$')

plt.show()
