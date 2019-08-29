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

dt = 0.01
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
    evolve = np.array([0., 0., 0., 0.])
    evolve[0] = d[0] + dt*d[2]
    evolve[1] = d[1] + dt*d[3]
    evolve[2] = d[2] + dt*(-d[0]*(1 + 2*v[1]) - 2*v[0]*d[1])
    evolve[3] = d[3] + dt*(-2*v[0]*d[0] + d[1]*(2*v[1] - 1))
    
    return evolve

H0 = 1/8
H = H0

xi = 0.0
#yi = 0.1    #Regular orbit
yi = -0.1   #Chaotic orbit
pyi = 0.0
pxi = sqrt(2*(H0 + p(yi, 3)/3 - p(xi, 2)*yi) - p(xi, 2) - p(yi, 2) - p(pyi, 2))

initConds = np.array([xi, yi, pxi, pyi])

i = 1
j = 0

np.random.seed(0)
rndm = np.array([np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None)])
dev1 = np.array([1.0, 0.0, 0.0, 0.0]) + rndm
dev2 = np.array([0.0, 1.0, 0.0, 0.0]) + rndm
dev3 = np.array([0.0, 0.0, 1.0, 0.0]) + rndm
dev4 = np.array([0.0, 0.0, 0.0, 1.0]) + rndm

mag1 = sqrt(p(dev1[0], 2) + p(dev1[1], 2) + p(dev1[2], 2) + p(dev1[3], 2))
mag2 = sqrt(p(dev2[0], 2) + p(dev2[1], 2) + p(dev2[2], 2) + p(dev2[3], 2))
mag3 = sqrt(p(dev3[0], 2) + p(dev3[1], 2) + p(dev3[2], 2) + p(dev3[3], 2))
mag4 = sqrt(p(dev4[0], 2) + p(dev4[1], 2) + p(dev4[2], 2) + p(dev4[3], 2))
dev1 = dev1/mag1
dev2 = dev2/mag2
dev3 = dev3/mag3
dev4 = dev4/mag4

print 'ICs: ', initConds
print 'Deviation vector 1: ', dev1
print 'Deviation vector 2: ', dev2
print 'Deviation vector 3: ', dev3
print 'Deviation vector 4: ', dev4
print

V = initConds

renorm = 100
rt = renorm*dt
k = 0

X = np.zeros((int(steps/renorm), 4))
X[0, :] = 1
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
        dev1 = evolveDev(V, dev1)
        dev2 = evolveDev(V, dev2)
        dev3 = evolveDev(V, dev3)
        dev4 = evolveDev(V, dev4)

        k1 = dt*func(V)
        k2 = dt*func(V + k1/2)
        k3 = dt*func(V + k2/2)
        k4 = dt*func(V + k3)
    
        V[0] = V[0] + (1.0/6.0 * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]))
        V[1] = V[1] + (1.0/6.0 * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]))
        V[2] = V[2] + (1.0/6.0 * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]))
        V[3] = V[3] + (1.0/6.0 * (k1[3] + 2*k2[3] + 2*k3[3] + k4[3]))
        
        k += 1

    mag1 = sqrt(p(dev1[0], 2) + p(dev1[1], 2) + p(dev1[2], 2) + p(dev1[3], 2))
    X[i, 0] = (i*rt*X[i - 1, 0]) + np.log(mag1)
    X[i, 0] = X[i, 0]/rt/(i + 1)
    dev1 = dev1/mag1

    dev2 = dev2 - np.dot(dev1, dev2)*dev1
    mag2 = sqrt(p(dev2[0], 2) + p(dev2[1], 2) + p(dev2[2], 2) + p(dev2[3], 2))
    X[i, 1] = (i*rt*X[i - 1, 1]) + np.log(mag2)
    X[i, 1] = X[i, 1]/rt/(i + 1)
    dev2 = dev2/mag2

    dev3 = dev3 - np.dot(dev1, dev3)*dev1 - np.dot(dev2, dev3)*dev2
    mag3 = sqrt(p(dev3[0], 2) + p(dev3[1], 2) + p(dev3[2], 2) + p(dev3[3], 2))
    X[i, 2] = (i*rt*X[i - 1, 2]) + np.log(mag3)
    X[i, 2] = X[i, 2]/rt/(i + 1)
    dev3 = dev3/mag3

    dev4 = dev4 - np.dot(dev1, dev4)*dev1 - np.dot(dev2, dev4)*dev2 - np.dot(dev3, dev4)*dev3
    mag4 = sqrt(p(dev4[0], 2) + p(dev4[1], 2) + p(dev4[2], 2) + p(dev4[3], 2))
    X[i, 3] = (i*rt*X[i - 1, 3]) + np.log(mag4)
    X[i, 3] = X[i, 3]/rt/(i + 1)
    dev4 = dev4/mag4
    
    i += 1

t1 = time.time()

print 'Final X1 estimate: ', X[int(steps/renorm) - 1, 0]
print 'Final X2 estimate: ', X[int(steps/renorm) - 1, 1]
print 'Final X3 estimate: ', X[int(steps/renorm) - 1, 2]
print 'Final X4 estimate: ', X[int(steps/renorm) - 1, 3] 
print
print
print 'Time elapsed: ', t1 - t0
print 'H = ', H0
if iterations > 0:
    print 'Iterations: ', iterations
else:
    print 'Iterations: ', steps

times = np.linspace(0, t_n, int(steps/renorm))

plt.plot(times, X[:, 0], 'bd--', linewidth = 0.5, markersize = 2, label = '$X_1$')
plt.plot(times, X[:, 1], 'rd--', linewidth = 0.5, markersize = 2, label = '$X_2$')
plt.plot(times, X[:, 2], 'gd--', linewidth = 0.5, markersize = 2, label = '$X_3$')
plt.plot(times, X[:, 3], 'kd--', linewidth = 0.5, markersize = 2, label = '$X_4$')
plt.xlabel('$t$')
plt.ylabel('$X_i$')
plt.xscale('log')
#plt.ylim(-0.2, 0.6)
plt.legend()

plt.show()
