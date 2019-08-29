from __future__ import division
from math import sqrt
from math import exp
from math import pow as p
import numpy as np
import matplotlib.pyplot as plt

V = np.zeros(2)
V_new = np.zeros(2)

K = 0.8#0.971635

tp = 2*np.pi

steps = 30

linIC_density = 20
initials = int(p(linIC_density, 2))

initConds = np.zeros([initials, 2])

for j in range(0, initials):
    initConds[j, 0] = (tp*j/linIC_density)%tp
    initConds[j, 1] = tp*np.floor(j/linIC_density)/linIC_density

def standMap(v):
    V_new[1] = v[1] + K*np.sin(v[0])
    V_new[0] = v[0] + V_new[1]
    
    return V_new

plt.scatter(initConds[:, 0], initConds[:, 1], s = 2, c = 'black')

for i in range(0, int(steps)):
    for k in range(0, initials):
        V = initConds[k, :]
        V = standMap(V)
        
        V[0] = V[0]%tp
        V[1] = V[1]%tp

        initConds[k, :] = V
        
        plt.scatter(V[0], V[1], s = 2, c = 'black')

plt.title('2D Standard Map with $K=0.971635$')
plt.xlabel('$x_1$')
plt.ylabel('$x_2$')

plt.show()
