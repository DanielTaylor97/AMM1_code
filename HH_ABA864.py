from __future__ import division
from math import sqrt
from math import exp
from math import pow as p
import numpy as np
import matplotlib.pyplot as plt
import time

rk_f = np.array([0.0, 0.0, 0.0, 0.0])
symp_int = np.array([0.0, 0.0, 0.0, 0.0])
interpVec = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
T = 0.0

dt = 0.1917
t_n = 5e4
V = np.array([0., 0., 0., 0.])
fn = V
steps = t_n/dt

H0 = 1/8
H = H0
initials = 1

def func(v):
    global rk_f
    rk_f[0] = v[2]                       #corresponds to t-deriv. of x
    rk_f[1] = v[3]                       #corresponds to t-deriv. of y
    rk_f[2] = -(v[0] + 2*v[0]*v[1])             #corresponds to t-deriv. of p_x
    rk_f[3] = v[1]*v[1] - v[1] - v[0]*v[0]      #corresponds to t-deriv. of p_y
    return rk_f
    # v = (x, y, px, py)

def e_LA(v, tau):
    symp_int[0] = v[0] + tau*v[2]
    symp_int[1] = v[1] + tau*v[3]
    symp_int[2] = v[2]
    symp_int[3] = v[3]

    return symp_int

def e_LB(v, tau):
    symp_int[0] = v[0]
    symp_int[1] = v[1]
    symp_int[2] = v[2] - tau*(v[0] + 2*v[0]*v[1])
    symp_int[3] = v[3] - tau*(v[1] + p(v[0], 2) - p(v[1], 2))

    return symp_int

def interp(vn, tn):
    global interpVec
    tn = tn*dt
    T = -vn[0]
    f = func(vn)
    interpVec[0] = 0
    interpVec[1] = vn[1] + T*(f[1])/(f[0])
    interpVec[2] = vn[2] + T*(f[2])/(f[0])
    interpVec[3] = vn[3] + T*(f[3])/(f[0])
    interpVec[4] = tn + T/f[0]
    return interpVec

t0 = time.time()

initConds = [[0., 0.]]
count = 0
condArray = [0]

f= open("initCondsContour.txt","r")
#f= open("initCondsH=0125.txt","r")

if f.mode == 'r':
    contents = f.readlines()
    for line in contents:
        count += 1
        condArray = np.append(condArray, float(line))

for s in range(0, int(count/2)):
    initConds = np.append(initConds, [[condArray[2*s + 1], condArray[2*s + 2]]], axis = 0)

initConds = np.delete(initConds, 0, 0)
#initConds = np.array([0.0, 0.1, 0.0, 0.0])   #Regular Orbit
#initConds = np.array([0.0, -0.1, 0.0, 0.0])  #Chaotic Orbit
initConds = np.array([0.0, 0.44, 0.0, 0.13])  #Sticky Orbit

#initials = np.size(initConds, 0)
initials = 1
print(initials)

valsX = np.zeros(initials)
valsY = np.array([initConds[1]])#initConds[:, 0]
valsPy = np.array([initConds[3]])#initConds[:, 1]
valsPx = np.zeros(initials)

xi = 0.0
yi = 0.0
pxi = 0.0
pyi = 0.0

for k in range(0, initials):
    xi = valsX[k]
    yi = valsY[k]
    pyi = valsPy[k]
    #print '', xi, ', ', yi, ', ', pyi
    pxi = sqrt(2*(H0 + p(yi, 3)/3 - p(xi, 2)*yi) - p(xi, 2) - p(yi, 2) - p(pyi, 2))
    valsPx[k] = pxi
    #print(pxi)

initConds[2] = pxi
#print(initConds)

plt.scatter(valsY, valsPy, s = 1, c = 'black')

i = 1
j = 0

a1 = 0.0711334264982231177779387300061549964174
a2 = 0.241153427956640098736487795326289649618
a3 = 0.521411761772814789212136078067994229991
a4 = -0.333698616227678005726562603400438876027
b1 = 0.183083687472197221961703757166430291072
b2 = 0.310782859898574869507522291054262796375
b3 = -0.0265646185119588006972121379164987592663
b4 = 0.0653961422823734184559721793911134363710

while i <= (steps - 1.0):
    #l=0
    for l in range(0, initials):
        V[0] = valsX[l]
        V[1] = valsY[l]
        V[2] = valsPx[l]
        V[3] = valsPy[l]
        
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

        if valsX[l]*V[0] <= 0:# and valsX[l] > V[0]:
            fn = interp(V, i)
            plt.scatter(fn[1], fn[3], s=1, c='black')
            j = j + 1

        valsX[l] = V[0]
        valsY[l] = V[1]
        valsPx[l] = V[2]
        valsPy[l] = V[3]
        
        #l = l + 1
    
    i = i + 1 

t1 = time.time()

print
print 'time: ', t1 - t0
print 'H = ', H0
print 'iterations: ', steps
print 'crossings: ', j

plt.xlabel('$y$', fontsize = 12)
plt.ylabel('$p_y$', fontsize = 12)

plt.show()
