from __future__ import division
from math import sqrt
from math import exp
from math import pow as p
import numpy as np
import matplotlib.pyplot as plt
import time

V = np.zeros(2)
V_new = np.zeros(2)

K = 0.971635

tp = 2*np.pi

steps = 10000

#linIC_density = 20
#initials = int(p(linIC_density, 2))
#
#initConds = np.zeros([initials, 2])
#
#for j in range(0, initials):
#    initConds[j, 0] = (j/linIC_density)%1.0
#    initConds[j, 1] = np.floor(j/linIC_density)/linIC_density

def standMap(v):
    V_new[1] = v[1] + K*np.sin(tp*v[0])/tp
    V_new[0] = v[0] + V_new[1]
    
    return V_new

#initConds = np.array([0.34, 0.42])  #Regular Orbit
initConds = np.array([0.03, 0.88])  #Chaotic Orbit
initials = 2

plt.scatter(initConds[0], initConds[1], s = 2, c = 'black')

#for i in range(0, int(steps)):
#    for k in range(0, initials):
#        V = initConds[k, :]
#        V = standMap(V)
#        
#        V[0] = V[0]%1.0
#        V[1] = V[1]%1.0
#
#        initConds[k, :] = V
#        
#        plt.scatter(V[0], V[1], s = 2, c = 'black')

t0 = time.time()

V = initConds
for i in range(0, int(steps)):
    V = standMap(V)

    V[0] = V[0]%1.0
    V[1] = V[1]%1.0

    plt.scatter(V[0], V[1], s = 2, c = 'black')

t1 = time.time()

print 'Time elapsed: ', t1 - t0

plt.xlabel('$x_1$')
plt.ylabel('$x_2$')

plt.show()
