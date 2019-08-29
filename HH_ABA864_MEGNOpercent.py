from __future__ import division
from math import sqrt
from math import exp
from math import pow as p
import numpy as np
import matplotlib.pyplot as plt
import time

##BEFORE RUNNING, COMMENT OUT TEXT FILE WRITING!!!

#0.5004 for 60*60 grid, t_n=10^3, thresh=-2.5 (1006s)

rk_f = np.array([0., 0., 0., 0., 0., 0., 0., 0.])
symp_int = np.array([0., 0., 0., 0., 0., 0., 0., 0.])
interpVec = np.array([0., 0., 0., 0., 0., 0., 0., 0.])
T = 0.0

dt = 0.1917
t_n = 1e3
V = np.array([0., 0., 0., 0., 0., 0., 0., 0.])
fn = V
steps = t_n/dt

m = 1
n = -1

H0 = 1/8
H = H0




initConds1 = [[0., 0.]]
count = 0
condArray = [0]

#f= open("initCondsH=0125.txt","r")
f= open("initCondsContour.txt","r")

if f.mode == 'r':
    contents = f.readlines()
    for line in contents:
        count += 1
        condArray = np.append(condArray, float(line))

for s in range(0, int(count/2)):
    initConds1 = np.append(initConds1, [[condArray[2*s + 1], condArray[2*s + 2]]], axis = 0)

linIC_density = 60  #200 will take 3hrs

initConds1 = np.delete(initConds1, 0, 0)
initials = np.size(initConds1, 0)
initConds = np.zeros([initials, 4])
print 'ICs: ', initials




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

i = 1

a1 = 0.0711334264982231177779387300061549964174
a2 = 0.241153427956640098736487795326289649618
a3 = 0.521411761772814789212136078067994229991
a4 = -0.333698616227678005726562603400438876027
b1 = 0.183083687472197221961703757166430291072
b2 = 0.310782859898574869507522291054262796375
b3 = -0.0265646185119588006972121379164987592663
b4 = 0.0653961422823734184559721793911134363710

renorm = 100
rt = renorm*dt
k = 0

Y = np.zeros(int(steps/renorm))
Y_ave = np.zeros(int(steps/renorm))
Y[0] = 1
Y_ave[0] = 1

vt = np.zeros(2)

for w in range(0, initials):
    vt = np.array([initConds1[w, 0], initConds1[w, 1]])
    initConds[w, 0] = 0
    initConds[w, 1] = vt[0]
    initConds[w, 3] = vt[1]
    initConds[w, 2] = sqrt(2*(H0 + p(vt[0], 3)/3) - p(vt[0], 2) - p(vt[1], 2))

print 'GO:'

t0 = time.time()

mGrid = np.zeros(initials)
times = np.linspace(0, steps, int(steps/renorm))

for l in range(0, initials):
    dev = np.array([np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None), np.random.uniform(low = -1.0, high = 1.0, size = None)])
    mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))
    dev = dev/mag
    mag = sqrt(p(dev[0], 2) + p(dev[1], 2) + p(dev[2], 2) + p(dev[3], 2))

    V = np.concatenate((initConds[l, :], dev))

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

    A = np.vstack([times, np.ones(len(times))]).T
    a, b = np.linalg.lstsq(A, 2*Y_ave, rcond = None)[0]
    mGrid[l] = a

    if l%10 == 0:
        print(l/initials)



t1 = time.time()

print
print 'Time elapsed: ', t1 - t0
print 'H = ', H0
print 'Iterations: ', steps



#f= open("HHcontour.txt","w+")
#
#for k in range(0, initials):
#    f.write(str(mGrid[k]))
#    f.write("\n")
#
#f.close()



grid = np.zeros((linIC_density, linIC_density))
chaotic = 0
total = initials#p(linIC_density, 2)
count = 0

Y = np.linspace(-0.4397, 0.6736, linIC_density)
PY = np.linspace(-0.5, 0.5, linIC_density)

for y in range(0, linIC_density):
    for py in range(0, linIC_density):
        if Y[y] == initConds[count, 1] and PY[py] == initConds[count, 3]:
            grid[py, y] = np.log10(np.abs(mGrid[count]))
            count += 1
            
            if grid[py, y] >-2.5:
                chaotic += 1
        else:
            grid[py, y] = np.nan
        if count == initials:
            break
    if count == initials:
        break

print
print 'Proportion of orbits that are chaotic: ', chaotic/total

plt.contourf(Y, PY, grid)
plt.colorbar(orientation='horizontal')
plt.xlabel('$y$')
plt.ylabel('$p_y$')

plt.show()
