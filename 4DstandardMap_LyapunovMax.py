from __future__ import division
from math import sqrt
from math import exp
from math import pow as p
import numpy as np
import matplotlib.pyplot as plt

initConds = np.array([0.55, 0.1, 0.54, 0.01])       #Regular Orbit

V_new = np.zeros(8)

K = 0.5
B = 0.05

c1 = 0
c3 = 0
c13 = 0

tp = 2*np.pi

steps = 1e6

def standMap(v):
    V_new = np.zeros(8)

    V_new[0] = v[0] + v[1] + K*(np.sin(tp*v[0]))/tp - B*(np.sin(tp*(v[2] - v[0])))/tp
    V_new[1] = V_new[0] - v[0]
    V_new[2] = v[2] + v[3] + K*(np.sin(tp*v[2]))/tp - B*(np.sin(tp*(v[0] - v[2])))/tp
    V_new[3] = V_new[2] - v[2]

    c1 = K*np.cos(tp*v[0])
    c3 = K*np.cos(tp*v[2])
    c13 = B*np.cos(tp*(v[0] - v[2]))

    V_new[5] = v[5] + v[4]*c1 + c13*(v[4] - v[6])
    V_new[4] = v[4] + V_new[5]
    V_new[7] = v[7] + v[6]*c3 + c13*(v[6] - v[4])
    V_new[6] = v[6] + V_new[7]

    return V_new

i = 1

#np.random.seed(0)
dev = np.array([np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None)])
mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))
dev = dev/mag
mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))

print 'Regular orbit ICs: ', initConds
print 'Regular orbit deviation vector: ', dev
print 'Magnitude: ', mag
print

V = np.concatenate((initConds, dev))

renorm = 10
k = 0

X1 = np.zeros(int(steps/renorm))
X1[0] = 1
times = np.linspace(0, steps, int(steps/renorm))

while i < steps/renorm:
    k = 0
    while k < renorm:
        V = standMap(V)
        
        V[0] = V[0]%1.0
        V[1] = V[1]%1.0
        V[2] = V[2]%1.0
        V[3] = V[3]%1.0

        k += 1

    dev = np.array([V[4], V[5], V[6], V[7]])
    mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))
    
    X1[i] = (i*renorm*X1[i - 1]) + np.log(mag)
    X1[i] = X1[i]/renorm/(i + 1)
    dev = dev/mag

    V[4] = dev[0]
    V[5] = dev[1]
    V[6] = dev[2]
    V[7] = dev[3]

    i += 1

times = np.linspace(0, steps, int(steps/renorm))
plt.plot(times, X1, 'bd--', linewidth = 0.5, markersize = 2, label = 'Regular Orbit')

print 'Regular orbit final mLE approximation: ', X1[int(steps/renorm) - 1]
print

initConds = np.array([0.55, 0.1, 0.005, 0.01])      #Chaotic Orbit

#np.random.seed(0)
dev = np.array([np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None)])
mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))
dev = dev/mag
mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))

print 'Chaotic orbit ICs: ', initConds
print 'Chaotic orbit deviation vector: ', dev
print 'Magnitude: ', mag
print

V = np.concatenate((initConds, dev))
i = 1

while i < steps/renorm:
    k = 0
    while k < renorm:
        V = standMap(V)
        
        V[0] = V[0]%1.0
        V[1] = V[1]%1.0
        V[2] = V[2]%1.0
        V[3] = V[3]%1.0

        k += 1

    dev = np.array([V[4], V[5], V[6], V[7]])
    mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))
    
    X1[i] = (i*renorm*X1[i - 1]) + np.log(mag)
    X1[i] = X1[i]/renorm/(i + 1)
    dev = dev/mag

    V[4] = dev[0]
    V[5] = dev[1]
    V[6] = dev[2]
    V[7] = dev[3]

    i += 1

print 'Chaotic orbit final mLE estimate: ', X1[int(steps/renorm) - 1]

plt.plot(times, X1, 'rd--', linewidth = 0.5, markersize = 2, label = 'Chaotic Orbit')
plt.xlabel('$t$')
plt.ylabel('$mLE$')
plt.xscale('log')
plt.yscale('log')
plt.legend()

plt.show()
