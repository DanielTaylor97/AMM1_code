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

steps = 1e4

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

X1 = np.zeros(int(steps/renorm))
X1[0] = 1

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
    
        X1[i] = (i*renorm*X1[i - 1]) + np.log(mag)
        X1[i] = X1[i]/(i + 1)/renorm
        dev = dev/mag
    
        V[2] = dev[0]
        V[3] = dev[1]
        
        i += 1

    mGrid[l] = X1[np.size(X1) - 1]

t1 = time.time()
print 'Time Elapsed: ', t1 - t0

grid = np.zeros((linIC_density, linIC_density))
chaotic = 0
total = p(linIC_density, 2)

for x in range(0, linIC_density):
    for y in range(0, linIC_density):
        grid[x, y] = np.log10(np.abs(mGrid[x*linIC_density + y]))

        if grid[x, y] >-2.2:
            chaotic += 1

print 'Proportion of orbits that are chaotic: ', chaotic/total

x1 = np.linspace(0, 2*np.pi, linIC_density)
x2 = np.linspace(0, 2*np.pi, linIC_density)

plt.contourf(x1, x2, grid)
plt.colorbar(orientation='horizontal')
plt.xlabel('$x_1$')
plt.ylabel('$x_2$')

plt.show()
