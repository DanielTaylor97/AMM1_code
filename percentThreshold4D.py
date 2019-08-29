from __future__ import division
from math import sqrt
from math import exp
from math import pow as p
import numpy as np
import matplotlib.pyplot as plt
import time

t0 = time.time()

chaotic = 0

numPoints = 500
threshold = np.linspace(-4.0, 0, numPoints)
prevalence = np.zeros(numPoints)
#print(threshold)



sigma = [[0.0]]
count = 0

f = open("4DSMcontour.txt", "r")

if f.mode == "r":
    contents = f.readlines()
    for line in contents:
        count += 1
        sigma = np.append(sigma, float(line))

sigma = np.delete(sigma, 0, 0)
print(count)



for i in range(0, numPoints):
    chaotic = 0
    for j in range(0, count):
        if sigma[j] > threshold[i]:
            chaotic += 1
    prevalence[i] = chaotic/count



t1 = time.time()
print 'Time Elapsed: ', t1 - t0

#print(prevalence)
plt.plot(threshold, prevalence)
plt.xlabel('Threshold Value')
plt.ylabel('Estimated Proportion of Chaotic Orbits')

plt.show()
