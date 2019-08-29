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

initConds = np.array([2.51257, 2.77352])    #Regular orbit

def standMap(v):
    V_new[1] = v[1] + K*np.sin(tp*v[0])/tp
    V_new[0] = v[0] + V_new[1]

    V_new[3] = v[2]*K*np.cos(tp*v[0])/tp + v[3]
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

renorm = 10
k = 0

Y = np.zeros(int(steps/renorm))
Y_ave = np.zeros(int(steps/renorm))
Y[0] = 1
Y_ave[0] = 1

t0 = time.time()
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

times = np.linspace(0, steps, int(steps/renorm))
plt.plot(times, Y, 'bd--', linewidth = 0.5, markersize = 2, label = 'Regular Orbit $Y$')
plt.plot(times, Y_ave, 'gd--', linewidth = 0.5, markersize = 2, label = 'Regular Orbit $<Y>$')




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

t1 = time.time()

print
print

print 'Time elapsed: ', t1 - t0
print 'Iterations: ', steps

plt.plot(times, Y, 'rd--', linewidth = 0.5, markersize = 2, label = 'Chaotic Orbit $Y$')
plt.plot(times, 2*Y_ave, 'cd--', linewidth = 0.5, markersize = 2, label = 'Chaotic Orbit $2<Y>$')

A = np.vstack([times, np.ones(len(times))]).T
a, b = np.linalg.lstsq(A, 2*Y_ave, rcond = None)[0]
rLine = a*times + b
a1 = "{:.4E}".format(a)
plt.plot(times, rLine, 'k--', linewidth = 1.0, markersize = 2, label = 'C.O. $2<Y>$ linear reg. m=%s' %a1)
print 'a = ', a
print 'b = ', b

plt.xlabel('$t$')
plt.ylabel('$MEGNO$ Indicator Value')
#plt.xscale('log')
#plt.yscale('log')
plt.legend()

plt.show()
