from __future__ import division
from math import sqrt
from math import exp
from math import pow as p
import numpy as np
import matplotlib.pyplot as plt
import time

V = np.zeros(4)
V_new = np.zeros(6)

K = 0.971635

tp = 2*np.pi

steps = 1e6
thresholdX = 1e-8

#initConds = np.array([2.51257, 2.77352])    #Regular orbit
initConds = np.array([0.2, 5.5])    #Chaotic orbit

def standMap(v):
    V_new = np.zeros(6)
    V_new[1] = v[1] + K*np.sin(v[0])
    V_new[0] = v[0] + V_new[1]

    V_new[3] = v[2]*K*np.cos(v[0]) + v[3]
    V_new[2] = V_new[3] + v[2]

    V_new[5] = v[4]*K*np.cos(v[0]) + v[5]
    V_new[4] = V_new[5] + v[4]
    
    return V_new

#np.random.seed(0)
rndm = np.array([np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None)])
dev1 = np.array([1.0, 0.0]) + rndm
dev2 = np.array([0.0, 1.0]) + rndm

mag1 = sqrt(p(dev1[0], 2) + p(dev1[1], 2))
mag2 = sqrt(p(dev2[0], 2) + p(dev2[1], 2))
dev1 = dev1/mag1
dev2 = dev2/mag2

print 'ICs: ', initConds
print 'Deviation vector 1: ', dev1
print 'Deviation vector 2: ', dev2
print

V = np.concatenate((initConds, dev1))
V = np.concatenate((V, dev2))

renorm = 10
k = 0

X = np.zeros((int(steps/renorm), 2))
X[0, :] = 1

t0 = time.time()
i = 1

while i < int(steps/renorm):
    k = 0
    while k < renorm:
        V = standMap(V)
    
        V[0] = V[0]%tp
        V[1] = V[1]%tp

        k += 1

    dev1 = np.array([V[2], V[3]])
    mag1 = sqrt(p(dev1[0], 2) + p(dev1[1], 2))
    dev2 = np.array([V[4], V[5]])

    X[i, 0] = (i*renorm*X[i - 1, 0]) + np.log(mag1)
    X[i, 0] = X[i, 0]/(i + 1)/renorm
    dev1 = dev1/mag1

    dev2 = dev2 - np.dot(dev1, dev2)*dev1
    mag2 = sqrt(p(dev2[0], 2) + p(dev2[1], 2))
    X[i, 1] = (i*renorm*X[i - 1, 1]) + np.log(mag2)
    X[i, 1] = X[i, 1]/(i + 1)/renorm
    dev2 = dev2/mag2

    V[2] = dev1[0]
    V[3] = dev1[1]
    V[4] = dev2[0]
    V[5] = dev2[1]
    
    i += 1

t1 = time.time()

print 'Final X1 estimate: ', X[int(steps/renorm) - 1, 0]
print 'Final X2 estimate: ', X[int(steps/renorm) - 1, 1]
print
print
print 'Time elapsed: ', t1 - t0
print 'Iterations: ', steps

times = np.linspace(0, steps, int(steps/renorm))

#plt.plot(times, X[:, 0], 'bd--', linewidth = 0.5, markersize = 2, label = 'Vec1')
#plt.plot(times, X[:, 1], 'rd--', linewidth = 0.5, markersize = 2, label = 'Vec2')
plt.plot(times, np.abs(X[:, 0] + X[:, 1]), 'bd--', linewidth = 0.5, markersize = 2, label = '$X_1 + X_2$')
plt.xlabel('$t$')
plt.ylabel('Sum of $X_i$')
plt.xscale('log')
plt.yscale('log')
plt.legend()

plt.show()
