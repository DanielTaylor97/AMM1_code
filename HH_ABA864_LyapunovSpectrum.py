from __future__ import division
from math import sqrt
from math import exp
from math import pow as p
import numpy as np
import matplotlib.pyplot as plt
import time

rk_f = np.array([0., 0., 0., 0., 0., 0., 0., 0.])
symp_int = np.zeros(20)
interpVec = np.array([0., 0., 0., 0., 0., 0., 0., 0.])
T = 0.0

dt = 0.1917
t_n = 100000
V = np.array([0., 0., 0., 0., 0., 0., 0., 0.])
fn = V
steps = t_n/dt
thresholdX = 1e-8

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

def e_LA(v, tau):
    symp_int[0] = v[0] + tau*v[2]
    symp_int[1] = v[1] + tau*v[3]
    symp_int[2] = v[2]
    symp_int[3] = v[3]

    symp_int[4] = v[4] + tau*v[6]
    symp_int[5] = v[5] + tau*v[7]
    symp_int[6] = v[6]
    symp_int[7] = v[7]

    symp_int[8] = v[8] + tau*v[10]
    symp_int[9] = v[9] + tau*v[11]
    symp_int[10] = v[10]
    symp_int[11] = v[11]

    symp_int[12] = v[12] + tau*v[14]
    symp_int[13] = v[13] + tau*v[15]
    symp_int[14] = v[14]
    symp_int[15] = v[15]
    
    symp_int[16] = v[16] + tau*v[18]
    symp_int[17] = v[17] + tau*v[19]
    symp_int[18] = v[18]
    symp_int[19] = v[19]

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

    symp_int[8] = v[8]
    symp_int[9] = v[9]
    symp_int[10] = v[10] - tau*(v[8]*(1 + 2*v[1]) + 2*v[0]*v[9])
    symp_int[11] = v[11] + tau*(-2*v[0]*v[8] + v[9]*(-1 + 2*v[1]))

    symp_int[12] = v[12]
    symp_int[13] = v[13]
    symp_int[14] = v[14] - tau*(v[12]*(1 + 2*v[1]) + 2*v[0]*v[13])
    symp_int[15] = v[15] + tau*(-2*v[0]*v[12] + v[13]*(-1 + 2*v[1]))

    symp_int[16] = v[16]
    symp_int[17] = v[17]
    symp_int[18] = v[18] - tau*(v[16]*(1 + 2*v[1]) + 2*v[0]*v[17])
    symp_int[19] = v[19] + tau*(-2*v[0]*v[16] + v[17]*(-1 + 2*v[1]))

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

xi = 0.0
#yi = 0.1    #Regular orbit
#yi = -0.1   #Chaotic orbit
#pyi = 0.44

yi = 0.44   #Sticky
pyi = 0.13

pxi = sqrt(2*(H0 + p(yi, 3)/3 - p(xi, 2)*yi) - p(xi, 2) - p(yi, 2) - p(pyi, 2))

initConds = np.array([xi, yi, pxi, pyi])

i = 1

a1 = 0.0711334264982231177779387300061549964174
a2 = 0.241153427956640098736487795326289649618
a3 = 0.521411761772814789212136078067994229991
a4 = -0.333698616227678005726562603400438876027
b1 = 0.183083687472197221961703757166430291072
b2 = 0.310782859898574869507522291054262796375
b3 = -0.0265646185119588006972121379164987592663
b4 = 0.0653961422823734184559721793911134363710

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
rt = renorm*dt
k = 0

X = np.zeros((int(steps/renorm), 4))
X[0, :] = 1

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

    dev1 = np.array([V[4], V[5], V[6], V[7]])
    mag1 = sqrt(p(dev1[0], 2) + p(dev1[1], 2) + p(dev1[2], 2) + p(dev1[3], 2))
    dev2 = np.array([V[8], V[9], V[10], V[11]])
    dev3 = np.array([V[12], V[13], V[14], V[15]])
    dev4 = np.array([V[16], V[17], V[18], V[19]])

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

#plt.plot(times, X[:, 0], 'bd--', linewidth = 0.5, markersize = 2, label = 'Vec1')
#plt.plot(times, X[:, 1], 'rd--', linewidth = 0.5, markersize = 2, label = 'Vec2')
#plt.plot(times, X[:, 2], 'gd--', linewidth = 0.5, markersize = 2, label = 'Vec3')
#plt.plot(times, X[:, 3], 'kd--', linewidth = 0.5, markersize = 2, label = 'Vec4')
plt.plot(times, np.abs(X[:, 0] + X[:, 3]), 'bd--', linewidth = 0.5, markersize = 2, label = '$X_1 + X_4$')
plt.plot(times, np.abs(X[:, 1] + X[:, 2]), 'rd--', linewidth = 0.5, markersize = 2, label = '$X_2 + X_3$')
plt.xlabel('$t$')
plt.ylabel('Sum of $X_i$')
plt.xscale('log')
plt.yscale('log')
plt.legend()

plt.show()
