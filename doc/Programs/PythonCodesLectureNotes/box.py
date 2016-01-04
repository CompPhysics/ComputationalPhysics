#!/usr/bin/env python
from  matplotlib import pyplot as plt
from math import exp
import numpy as np
import random

# initial number of particles
N0 = 1000
MaxTime = 10*N0
values = np.zeros(MaxTime)   
time = np.zeros(MaxTime)   
random.seed() 
# initial number of particles in left half
nleft = N0
for t in range (0, MaxTime, 1):
    if N0*random.random() <= nleft: 
       nleft -= 1
    else: 
       nleft += 1
    time[t] = t
    values[t] = nleft

# Finally we plot the results
plt.plot(time, values,'b-')
plt.axis([0,MaxTime, N0/4, N0])
plt.xlabel('$t$')
plt.ylabel('$N$')
plt.title('Number of particles in left half')
plt.savefig('box.pdf')
plt.show()


