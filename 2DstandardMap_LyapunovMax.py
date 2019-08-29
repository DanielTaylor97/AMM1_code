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

steps = 1e6
thresholdX = 1e-8

initConds = np.array([2.51257, 2.77352])    #Regular orbit

def standMap(v):
    V_new[1] = v[1] + K*np.sin(v[0])
    V_new[0] = v[0] + V_new[1]

    V_new[3] = v[2]*K*np.cos(v[0]) + v[3]
    V_new[2] = V_new[3] + v[2]
    
    return V_new

#np.random.seed(0)
dev = np.array([np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None)])
mag = sqrt(p(dev[0], 2) + p(dev[1], 2))
dev = dev/mag
mag = sqrt(p(dev[0], 2) + p(dev[1], 2))

print 'Regular orbit ICs: ', initConds
print 'Regular orbit deviation vector: ', dev
print 'Magnitude: ', mag
print

V = np.concatenate((initConds, dev))

renorm = 100
k = 0

X1 = np.zeros(int(steps/renorm))
X1[0] = 1

t0 = time.time()
iterations = 0
i = 1

while i < int(steps/renorm):
    k = 0
    while k < renorm:
        V = standMap(V)
    
        V[0] = V[0]%tp
        V[1] = V[1]%tp

        k += 1

    dev = np.array([V[2], V[3]])
    mag = sqrt(p(dev[0], 2) + p(dev[1], 2))

    X1[i] = (i*renorm*X1[i - 1]) + np.log(mag)
    X1[i] = X1[i]/(i + 1)/renorm
    dev = dev/mag

#    if X1[i] < thresholdX:
#        print '!!! ', X1[i], ' less than ', thresholdX
#        iterations = i
#        i = int(steps)

    V[2] = dev[0]
    V[3] = dev[1]
    
    i += 1

times = np.linspace(0, steps, int(steps/renorm))
plt.plot(times, X1, 'bd--', linewidth = 0.5, markersize = 2, label = 'Regular Orbit')

print 'Regular orbit final mLE approximation: ', X1[int(steps/renorm) - 1]
print

initConds = np.array([0.2, 5.5])    #Chaotic orbit

#np.random.seed(0)
dev = np.array([np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None)])
mag = sqrt(p(dev[0], 2) + p(dev[1], 2))
dev = dev/mag
mag = sqrt(p(dev[0], 2) + p(dev[1], 2))

print 'Chaotic orbit ICs: ', initConds
print 'Chaotic orbit deviation vector: ', dev
print 'Magnitude: ', mag
print

V = np.concatenate((initConds, dev))
i = 1

while i < int(steps/renorm):
    k = 0
    while k < renorm:
        V = standMap(V)
    
        V[0] = V[0]%tp
        V[1] = V[1]%tp

        k += 1

    dev[0] = V[2]
    dev[1] = V[3]

    mag = sqrt(p(dev[0], 2) + p(dev[1], 2))

    X1[i] = (i*renorm*X1[i - 1]) + np.log(mag)
    X1[i] = X1[i]/(i + 1)/renorm
    dev = dev/mag

    if X1[i] < thresholdX:
        print '!!! ', X1[i], ' less than ', thresholdX
        iterations = i
        i = int(steps)

    V[2] = dev[0]
    V[3] = dev[1]
    i += 1

t1 = time.time()

print 'Chaotic orbit final mLE approximation: ', X1[int(steps/renorm) - 1]
print
print

print 'Time elapsed: ', t1 - t0
if iterations > 0:
    print 'Iterations: ', iterations
else:
    print 'Iterations: ', steps

plt.plot(times, X1, 'rd--', linewidth = 0.5, markersize = 2, label = 'Chaotic Orbit')
plt.xlabel('$t$')
plt.ylabel('$mLE$')
plt.xscale('log')
plt.yscale('log')
plt.legend()

plt.show()
