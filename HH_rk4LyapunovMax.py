from __future__ import division
from math import sqrt
from math import exp
from math import pow as p
import numpy as np
import matplotlib.pyplot as plt
import time

rk_f = np.array([0.0, 0.0, 0.0, 0.0])
evolve = np.array([0., 0., 0., 0.])
T = 0.0

dt = 0.001      ##NECESSARY FOR FINAL VALUES NOT TO TAIL OFF WEIRDLY
t_n = 10000
V = np.array([0., 0., 0., 0.])
fn = V
steps = t_n/dt

def func(v):
    global rk_f
    rk_f[0] = v[2]                       #corresponds to t-deriv. of x
    rk_f[1] = v[3]                       #corresponds to t-deriv. of y
    rk_f[2] = -(v[0] + 2*v[0]*v[1])             #corresponds to t-deriv. of p_x
    rk_f[3] = p(v[1], 2) - v[1] - p(v[0], 2)     #corresponds to t-deriv. of p_y
    return rk_f

def evolveDev(v, d):
    evolve[0] = d[0] + dt*d[2]
    evolve[1] = d[1] + dt*d[3]
    evolve[2] = d[2] + dt*(-d[0]*(1 + 2*v[1]) - 2*v[0]*d[1])
    evolve[3] = d[3] + dt*(-2*v[0]*d[0] + d[1]*(2*v[1] - 1))
    
    return evolve

H0 = 1/8
H = H0

xi = 0.0
yi = 0.1    #Regular orbit
pyi = 0.0
pxi = sqrt(2*(H0 + p(yi, 3)/3 - p(xi, 2)*yi) - p(xi, 2) - p(yi, 2) - p(pyi, 2))

initConds = np.array([xi, yi, pxi, pyi])

i = 1
j = 0

#np.random.seed(0)
dev = np.array([np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None)])
mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))
dev = dev/mag
mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))

print 'Regular orbit ICs: ', initConds
print 'Regular orbit deviation vector: ', dev
print 'Magnitude: ', mag
print

V = initConds

renorm = 100
rt = renorm*dt
k = 0

X1 = np.zeros(int(steps/renorm))
X1[0] = 1
times = np.linspace(0, t_n, int(steps/renorm))

t0 = time.time()
iterations = 0

k1 = np.zeros(4)
k2 = np.zeros(4)
k3 = np.zeros(4)
k4 = np.zeros(4)

while i <= (steps/renorm - 1.0):
    k = 0
    while k < renorm:
        
        k1 = dt*func(V)
        k2 = dt*func(V + k1/2)
        k3 = dt*func(V + k2/2)
        k4 = dt*func(V + k3)

        dev = evolveDev(V, dev)
    
        V[0] = V[0] + (1.0/6.0 * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]))
        V[1] = V[1] + (1.0/6.0 * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]))
        V[2] = V[2] + (1.0/6.0 * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]))
        V[3] = V[3] + (1.0/6.0 * (k1[3] + 2*k2[3] + 2*k3[3] + k4[3]))
        
        k += 1
    
    mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))

    X1[i] = (i*rt*X1[i - 1]) + np.log(mag)
    X1[i] = X1[i]/rt/(i + 1)
    dev = dev/mag
    
    i += 1

times = np.linspace(0, t_n, int(steps/renorm))
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

V = initConds
i = 1

while i <= (steps/renorm - 1.0):
    k = 0
    while k < renorm:
        
        k1 = dt*func(V)
        k2 = dt*func(V + k1/2)
        k3 = dt*func(V + k2/2)
        k4 = dt*func(V + k3)

        dev = evolveDev(V, dev)
    
        V[0] = V[0] + (1.0/6.0 * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]))
        V[1] = V[1] + (1.0/6.0 * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]))
        V[2] = V[2] + (1.0/6.0 * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]))
        V[3] = V[3] + (1.0/6.0 * (k1[3] + 2*k2[3] + 2*k3[3] + k4[3]))
        
        k += 1
    
    mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))

    X1[i] = (i*rt*X1[i - 1]) + np.log(mag)
    X1[i] = X1[i]/rt/(i + 1)
    dev = dev/mag
    
    i += 1

t1 = time.time()

print 'Chaotic orbit final mLE approximation: ', X1[int(steps/renorm) - 1]
print
print
print 'Time elapsed: ', t1 - t0
print 'H = ', H0
print 'iterations: ', steps

plt.plot(times, X1, 'rd--', linewidth = 0.5, markersize = 2, label = 'Chaotic Orbit')
plt.xlabel('$t$')
plt.ylabel('$mLE$')
plt.xscale('log')
plt.yscale('log')
plt.legend()

plt.show()
