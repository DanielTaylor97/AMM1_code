from __future__ import division
from math import sqrt
from math import exp
from math import pow as p
import numpy as np
import matplotlib.pyplot as plt

H = 1/8

numPoints = 10000
y = np.linspace(-0.4396, 0.6735, numPoints)##
py = np.zeros(numPoints)

for i in range(0, numPoints):
    py[i] = sqrt(2*(H + p(y[i], 3)/3) - p(y[i], 2))##

initConds = [[0., 0.]]
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
    initConds = np.append(initConds, [[condArray[2*s + 1], condArray[2*s + 2]]], axis = 0)

initConds = np.delete(initConds, 0, 0)

#print(initConds)
print 'Initial Conditions: ', count/2

plt.scatter(initConds[:, 0], initConds[:, 1], s=1, c='black')
plt.plot(y, py, c='black')
plt.plot(y, -py, c='black')
plt.xlabel('y')
plt.ylabel('$p_y$')
plt.show()
