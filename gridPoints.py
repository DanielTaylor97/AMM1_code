from __future__ import division
from math import sqrt
from math import exp
from math import pow as p
import numpy as np
import matplotlib.pyplot as plt

H0 = 1/8
H = 0

length = 60     #Linear density of points

y = np.linspace(-0.4397, 0.6736, length)
py = np.linspace(-0.5, 0.5, length)

initConds = [[0., 0.]]

for i in range(0, length):#y
    for j in range(0, length):#py
        H = 0.5*(p(py[j], 2) + p(y[i], 2)) - p(y[i], 3)/3
        if H <= H0:
            initConds = np.append(initConds, [[y[i], py[j]]], axis = 0)

initConds = np.delete(initConds, 0, 0)
number = np.size(initConds, 0)
print 'ICs: ', number

#f= open("initCondsH=0125.txt","w+")
f= open("initCondsContour.txt","w+")

for k in range(0, number):
    for l in range(0, 2):
        f.write(str(initConds[k, l]))
        f.write("\n")

f.close()

#print(initConds)
#
#plt.scatter(initConds[:, 0], initConds[:, 1])
#plt.show()
