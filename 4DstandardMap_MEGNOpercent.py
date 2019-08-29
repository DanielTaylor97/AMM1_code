from __future__ import division
from math import sqrt
from math import exp
from math import pow as p
import numpy as np
import matplotlib.pyplot as plt
import time

##BEFORE RUNNING, COMMENT OUT TEXT FILE WRITING!!!

V_new = np.zeros(8)

K = 0.5
B = 0.05

c1 = 0
c3 = 0
c13 = 0

tp = 2*np.pi

steps = 1e3

m = 1
n = -1

linIC_density = 10
initials = int(p(linIC_density, 4))

initConds = np.zeros([initials, 4])

for j in range(0, initials):
    initConds[j, 0] = (j/linIC_density)%1.0
    initConds[j, 1] = (np.floor(j/linIC_density)/linIC_density)%1.0
    initConds[j, 2] = (np.floor(j/linIC_density/linIC_density)/linIC_density)%1.0
    initConds[j, 3] = (np.floor(j/linIC_density/linIC_density/linIC_density)/linIC_density)%1.0

#print(initConds)

def standMap(v):
    V_new = np.zeros(8)

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

    return V_new

i = 1

renorm = 10
k = 0

Y = np.zeros(int(steps/renorm))
Y_ave = np.zeros(int(steps/renorm))
Y[0] = 1
Y_ave[0] = 1

times = np.linspace(0, steps, int(steps/renorm))
A = np.vstack([times, np.ones(len(times))]).T
chaotic = 0
mGrid = np.zeros(initials)

t0 = time.time()




for l in range(0, initials):
    dev = np.array([np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None)])
    mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))
    dev = dev/mag
    mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))
    
    V = np.concatenate((initConds[l, :], dev))

    i = 1
    while i < steps/renorm:
        k = 0
        while k < renorm:
            V = standMap(V)
            
            V[0] = V[0]%1.0
            V[1] = V[1]%1.0
            V[2] = V[2]%1.0
            V[3] = V[3]%1.0
    
            k += 1
    
        dev = np.array([V[4], V[5], V[6], V[7]])
        mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))
        
        Y[i] = Y[i - 1]/(m + 1)/p(i, n)
        Y[i] = Y[i] + np.log(mag)*p(i + 1, m)
        Y[i] = (m + 1)*p(i + 1, n)*Y[i]
        dev = dev/mag
    
        Y_ave[i] = p(i, m + n + 1)*Y_ave[i - 1]
        Y_ave[i] = Y_ave[i] + Y[i]
        Y_ave[i] = Y_ave[i]/p(i + 1, m + n + 1)
    
        V[4] = dev[0]
        V[5] = dev[1]
        V[6] = dev[2]
        V[7] = dev[3]
    
        i += 1

    a, b = np.linalg.lstsq(A, 2*Y_ave, rcond = None)[0]
    a = np.log10(np.abs(a))

    mGrid[l] = a
    if a > -2.2:
        chaotic += 1
        
    #print((l + 1)/initials)
    




print 'Proportion of chaotic orbits: ', chaotic/initials
t1 = time.time()
print 'Time elapsed: ', t1 - t0



#f= open("4DSMcontour.txt","w+")
#
#for k in range(0, initials):
#    f.write(str(mGrid[k]))
#    f.write("\n")
#
#f.close()
