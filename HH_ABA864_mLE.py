from __future__ import division
from math import sqrt
from math import exp
from math import pow as p
import numpy as np
import matplotlib.pyplot as plt
import time

rk_f = np.array([0., 0., 0., 0., 0., 0., 0., 0.])
symp_int = np.array([0., 0., 0., 0., 0., 0., 0., 0.])
interpVec = np.array([0., 0., 0., 0., 0., 0., 0., 0.])
T = 0.0

dt = 0.1917
t_n = 1e5
V = np.array([0., 0., 0., 0., 0., 0., 0., 0.])
fn = V
steps = t_n/dt
thresholdX = 1e-8

H0 = 1/8
H = H0
initials = 1

def e_LA(v, tau):
    symp_int[0] = v[0] + tau*v[2]
    symp_int[1] = v[1] + tau*v[3]
    symp_int[2] = v[2]
    symp_int[3] = v[3]

    symp_int[4] = v[4] + tau*v[6]
    symp_int[5] = v[5] + tau*v[7]
    symp_int[6] = v[6]
    symp_int[7] = v[7]

    return symp_int

def e_LB(v, tau):
    symp_int[0] = v[0]
    symp_int[1] = v[1]
    symp_int[2] = v[2] - tau*(v[0] + 2*v[0]*v[1])
    symp_int[3] = v[3] - tau*(v[1] + p(v[0], 2) - p(v[1], 2))

    symp_int[4] = v[4]
    symp_int[5] = v[5]
    symp_int[6] = v[6] - tau*(v[4]*(1 + 2*v[1]) + 2*v[0]*v[5])
    symp_int[7] = v[7] + tau*(-2*v[0]*v[4] + v[5]*(-1 + 2*v[1]))

    return symp_int

xi = 0.0
yi = 0.1    #Regular orbit
pyi = 0.0
pxi = sqrt(2*(H0 + p(yi, 3)/3 - p(xi, 2)*yi) - p(xi, 2) - p(yi, 2) - p(pyi, 2))

initConds = np.array([xi, yi, pxi, pyi])

i = 1
j = 1

a1 = 0.0711334264982231177779387300061549964174
a2 = 0.241153427956640098736487795326289649618
a3 = 0.521411761772814789212136078067994229991
a4 = -0.333698616227678005726562603400438876027
b1 = 0.183083687472197221961703757166430291072
b2 = 0.310782859898574869507522291054262796375
b3 = -0.0265646185119588006972121379164987592663
b4 = 0.0653961422823734184559721793911134363710

#np.random.seed(0)
#dev = np.array([np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None)])
dev = np.array([1.0, 0.0, 0.0, 0.0])
mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))
dev = dev/mag
mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))

print 'Regular orbit ICs: ', initConds
print 'Regular orbit deviation vector: ', dev
print 'Magnitude: ', mag
print

V = np.concatenate((initConds, dev))

renorm = 100
rt = renorm*dt
k = 0

X1 = np.zeros(int(steps/renorm))
X1[0] = 1

t0 = time.time()
iterations = 0

while i <= (steps/renorm - 1.0):
    k = 0
    while k < renorm:
        V = e_LA(V, a1*dt)
        V = e_LB(V, b1*dt)
        V = e_LA(V, a2*dt)
        V = e_LB(V, b2*dt)
        V = e_LA(V, a3*dt)
        V = e_LB(V, b3*dt)
        V = e_LA(V, a4*dt)
        V = e_LB(V, b4*dt)
        V = e_LA(V, a4*dt)
        V = e_LB(V, b3*dt)
        V = e_LA(V, a3*dt)
        V = e_LB(V, b2*dt)
        V = e_LA(V, a2*dt)
        V = e_LB(V, b1*dt)
        V = e_LA(V, a1*dt)

        k += 1

    dev = np.array([V[4], V[5], V[6], V[7]])
    mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))

    X1[i] = (i*rt*X1[i - 1]) + np.log(mag)
    X1[i] = X1[i]/rt/(i + 1)
    dev = dev/mag

    if X1[i] < thresholdX:
        print '!!! ', X1[i], ' less than ', thresholdX
        iterations = i
        i = int(steps)

    V[4] = dev[0]
    V[5] = dev[1]
    V[6] = dev[2]
    V[7] = dev[3]
    
    i += 1

times = np.linspace(0, t_n, int(steps/renorm))
subtimes = times[int(1e3/renorm):]
plt.plot(times, X1, 'bd--', linewidth = 0.5, markersize = 2, label = 'Regular Orbit')

print 'Regular orbit final mLE approximation: ', X1[int(steps/renorm) - 1]
print

yi = -0.1   #Chaotic orbit
initConds = np.array([xi, yi, pxi, pyi])

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

while i <= (steps/renorm - 1.0):
    k = 0
    while k < renorm:
        V = e_LA(V, a1*dt)
        V = e_LB(V, b1*dt)
        V = e_LA(V, a2*dt)
        V = e_LB(V, b2*dt)
        V = e_LA(V, a3*dt)
        V = e_LB(V, b3*dt)
        V = e_LA(V, a4*dt)
        V = e_LB(V, b4*dt)
        V = e_LA(V, a4*dt)
        V = e_LB(V, b3*dt)
        V = e_LA(V, a3*dt)
        V = e_LB(V, b2*dt)
        V = e_LA(V, a2*dt)
        V = e_LB(V, b1*dt)
        V = e_LA(V, a1*dt)

        k += 1

    dev = np.array([V[4], V[5], V[6], V[7]])
    mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))

    X1[i] = (i*rt*X1[i - 1]) + np.log(mag)
    X1[i] = X1[i]/rt/(i + 1)
    dev = dev/mag

    if X1[i] < thresholdX:
        print '!!! ', X1[i], ' less than ', thresholdX
        iterations = i
        i = int(steps)

    V[4] = dev[0]
    V[5] = dev[1]
    V[6] = dev[2]
    V[7] = dev[3]
    
    i += 1

times = np.linspace(0, t_n, int(steps/renorm))
plt.plot(times, X1, 'rd--', linewidth = 0.5, markersize = 2, label = 'Chaotic Orbit')

subX1 = X1[int(1e3/renorm):]
b = np.mean(subX1)
rLine = b*np.ones(np.size(subX1))
b_ = "{:.4E}".format(b)
#plt.plot(subtimes, rLine, 'k--', linewidth = 0.5, markersize = 2, label = 'Chaotic Orbit $\chi_1$ ave. = %s' %b_)
print 'b = ', b

print 'Chaotic orbit final mLE approximation: ', X1[int(steps/renorm) - 1]
print

xi = 0.0
yi = 0.44    #Sticky orbit
pyi = 0.13
pxi = sqrt(2*(H0 + p(yi, 3)/3 - p(xi, 2)*yi) - p(xi, 2) - p(yi, 2) - p(pyi, 2))
initConds = np.array([xi, yi, pxi, pyi])

#np.random.seed(0)
dev = np.array([np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None)])
mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))
dev = dev/mag
mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))

print 'Sticky orbit ICs: ', initConds
print 'Sticky orbit deviation vector: ', dev
print 'Magnitude: ', mag
print

V = np.concatenate((initConds, dev))
i = 1

while i <= (steps/renorm - 1.0):
    k = 0
    while k < renorm:
        V = e_LA(V, a1*dt)
        V = e_LB(V, b1*dt)
        V = e_LA(V, a2*dt)
        V = e_LB(V, b2*dt)
        V = e_LA(V, a3*dt)
        V = e_LB(V, b3*dt)
        V = e_LA(V, a4*dt)
        V = e_LB(V, b4*dt)
        V = e_LA(V, a4*dt)
        V = e_LB(V, b3*dt)
        V = e_LA(V, a3*dt)
        V = e_LB(V, b2*dt)
        V = e_LA(V, a2*dt)
        V = e_LB(V, b1*dt)
        V = e_LA(V, a1*dt)

        k += 1

    dev = np.array([V[4], V[5], V[6], V[7]])
    mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))

    X1[i] = (i*rt*X1[i - 1]) + np.log(mag)
    X1[i] = X1[i]/rt/(i + 1)
    dev = dev/mag

    if X1[i] < thresholdX:
        print '!!! ', X1[i], ' less than ', thresholdX
        iterations = i
        i = int(steps)

    V[4] = dev[0]
    V[5] = dev[1]
    V[6] = dev[2]
    V[7] = dev[3]
    
    i += 1



t1 = time.time()

print 'Sticky orbit final mLE approximation: ', X1[int(steps/renorm) - 1]
print
print
print 'Time elapsed: ', t1 - t0
print 'H = ', H0
if iterations > 0:
    print 'Iterations: ', iterations
else:
    print 'Iterations: ', steps

#plt.plot(times, X1, 'gd--', linewidth = 0.5, markersize = 2, label = 'Sticky Orbit')

Yval_linsp = np.linspace(2e-4, 1e-1, int(steps/renorm))
plotVec = np.ones(int(steps/renorm))
plt.plot(times, (3.16e-3)*plotVec, 'k--', linewidth = 0.5, markersize = 2)
plt.plot((1e4)*plotVec, Yval_linsp, 'k--', linewidth = 0.5, markersize = 2)

plt.xlabel('$t$')
plt.ylabel('$mLE$')
plt.xscale('log')
plt.yscale('log')
plt.legend()

plt.show()
