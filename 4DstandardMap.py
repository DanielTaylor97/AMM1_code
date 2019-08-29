from __future__ import division
from math import sqrt
from math import exp
from math import pow as p
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

#V = np.array([0.55, 0.1, 0.54, 0.01])   #Regular Orbit
#V = np.array([0.55, 0.1, 0.005, 0.01])  #Chaotic Orbit
V = np.array([0.1, 0.10001, 0.1, 0.1])      #Chaotic Orbit
V_new = np.zeros(4)

K = 0.5
B = 0.05

tp = 2*np.pi
V = V*tp/12

steps = 1e5#9.2e5

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def standMap(v):
    V_new = np.zeros(4)
    V_new[1] = v[1] + K*np.sin(tp*v[0])/tp - B*np.sin(tp*(v[2] - v[0]))/tp
    V_new[0] = v[0] + V_new[1]
    V_new[3] = v[3] + K*np.sin(tp*v[2])/tp - B*np.sin(tp*(v[0] - v[2]))/tp
    V_new[2] = v[2] + V_new[3]

    return V_new

t0 = time.time()

for i in range(0, int(steps)):
    V = standMap(V)

    V[0] = V[0]%1.0
    V[1] = V[1]%1.0
    V[2] = V[2]%1.0
    V[3] = V[3]%1.0
    colour = np.sin(tp*V[3])/2 + 0.5

    ax.scatter(V[0], V[1], V[2], s=2, c=(colour, 0.0, 0.5), depthshade=True)

#    if i > 7.64e5:
#        ax.scatter(V[0], V[1], V[2], s=2, c=(colour, 0.0, 0.5), depthshade=True)
    #plt.scatter(V[0], V[1], s=2, c='black')

t1 = time.time()

print 'Elapsed time: ', t1 - t0

#ax.set_title('4D Standard Map with $K = 0.5$, $B = 0.05$')
ax.set_xlabel('$x_1$')
ax.set_ylabel('$x_2$')
ax.set_zlabel('$x_3$')

plt.show()
