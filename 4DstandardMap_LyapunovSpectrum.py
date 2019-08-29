from __future__ import division
from math import sqrt
from math import exp
from math import pow as p
import numpy as np
import matplotlib.pyplot as plt

initConds = np.array([0.55, 0.1, 0.54, 0.01])       #Regular Orbit
#initConds = np.array([0.55, 0.1, 0.005, 0.01])      #Chaotic Orbit
V_new = np.zeros(20)

K = 0.5
B = 0.05

c1 = 0
c3 = 0
c13 = 0

tp = 2*np.pi

steps = 1e6

def standMap(v):
    V_new = np.zeros(20)

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

    V_new[9] = v[9] + v[8]*c1 + c13*(v[8] - v[10])
    V_new[8] = v[8] + V_new[9]
    V_new[11] = v[11] + v[10]*c3 + c13*(v[10] - v[8])
    V_new[10] = v[10] + V_new[11]

    V_new[13] = v[13] + v[12]*c1 + c13*(v[12] - v[14])
    V_new[12] = v[12] + V_new[13]
    V_new[15] = v[15] + v[14]*c3 + c13*(v[14] - v[12])
    V_new[14] = v[14] + V_new[15]

    V_new[17] = v[17] + v[16]*c1 + c13*(v[16] - v[18])
    V_new[16] = v[16] + V_new[17]
    V_new[19] = v[19] + v[18]*c3 + c13*(v[18] - v[16])
    V_new[18] = v[18] + V_new[19]

    return V_new

i = 1

#np.random.seed(0)
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

V = np.concatenate((initConds, dev1))
V = np.concatenate((V, dev2))
V = np.concatenate((V, dev3))
V = np.concatenate((V, dev4))

renorm = 10
k = 0

X = np.zeros((int(steps/renorm), 4))
X[0, :] = 1

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

    dev1 = np.array([V[4], V[5], V[6], V[7]])
    mag1 = sqrt(p(dev1[0], 2) + p(dev1[1], 2) + p(dev1[2], 2) + p(dev1[3], 2))
    dev2 = np.array([V[8], V[9], V[10], V[11]])
    dev3 = np.array([V[12], V[13], V[14], V[15]])
    dev4 = np.array([V[16], V[17], V[18], V[19]])

    X[i, 0] = (i*renorm*X[i - 1, 0]) + np.log(mag1)
    X[i, 0] = X[i, 0]/renorm/(i + 1)
    dev1 = dev1/mag1

    dev2 = dev2 - np.dot(dev1, dev2)*dev1
    mag2 = sqrt(p(dev2[0], 2) + p(dev2[1], 2) + p(dev2[2], 2) + p(dev2[3], 2))
    X[i, 1] = (i*renorm*X[i - 1, 1]) + np.log(mag2)
    X[i, 1] = X[i, 1]/renorm/(i + 1)
    dev2 = dev2/mag2

    dev3 = dev3 - np.dot(dev1, dev3)*dev1 - np.dot(dev2, dev3)*dev2
    mag3 = sqrt(p(dev3[0], 2) + p(dev3[1], 2) + p(dev3[2], 2) + p(dev3[3], 2))
    X[i, 2] = (i*renorm*X[i - 1, 2]) + np.log(mag3)
    X[i, 2] = X[i, 2]/renorm/(i + 1)
    dev3 = dev3/mag3

    dev4 = dev4 - np.dot(dev1, dev4)*dev1 - np.dot(dev2, dev4)*dev2 - np.dot(dev3, dev4)*dev3
    mag4 = sqrt(p(dev4[0], 2) + p(dev4[1], 2) + p(dev4[2], 2) + p(dev4[3], 2))
    X[i, 3] = (i*renorm*X[i - 1, 3]) + np.log(mag4)
    X[i, 3] = X[i, 3]/renorm/(i + 1)
    dev4 = dev4/mag4

    V[4] = dev1[0]
    V[5] = dev1[1]
    V[6] = dev1[2]
    V[7] = dev1[3]
    V[8] = dev2[0]
    V[9] = dev2[1]
    V[10] = dev2[2]
    V[11] = dev2[3]
    V[12] = dev3[0]
    V[13] = dev3[1]
    V[14] = dev3[2]
    V[15] = dev3[3]
    V[16] = dev4[0]
    V[17] = dev4[1]
    V[18] = dev4[2]
    V[19] = dev4[3]

    i += 1

print 'Final X1 estimate: ', X[int(steps/renorm) - 1, 0]
print 'Final X2 estimate: ', X[int(steps/renorm) - 1, 1]
print 'Final X3 estimate: ', X[int(steps/renorm) - 1, 2]
print 'Final X4 estimate: ', X[int(steps/renorm) - 1, 3] 
print
print 'Iterations: ', steps

times = np.linspace(0, steps, int(steps/renorm))

#plt.plot(times, X[:, 0], 'bd--', linewidth = 0.5, markersize = 2, label = '$X_1$')
#plt.plot(times, X[:, 1], 'rd--', linewidth = 0.5, markersize = 2, label = '$X_2$')
#plt.plot(times, X[:, 2], 'gd--', linewidth = 0.5, markersize = 2, label = '$X_3$')
#plt.plot(times, X[:, 3], 'kd--', linewidth = 0.5, markersize = 2, label = '$X_4$')
plt.plot(times, np.abs(X[:, 0] + X[:, 3]), 'bd--', linewidth = 0.5, markersize = 2, label = '$X_1 + X_4$')
plt.plot(times, np.abs(X[:, 1] + X[:, 2]), 'rd--', linewidth = 0.5, markersize = 2, label = '$X_2 + X_3$')
plt.xlabel('$t$')
plt.ylabel('Sum of $X_i$')
plt.xscale('log')
plt.yscale('log')
plt.legend()

plt.show()
