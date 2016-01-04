#!/usr/bin/env python
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
import random

# number of random walkers
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


#!/usr/bin/env python
import numpy as np

import random

# initialize the rng with a seed
random.seed() 
counts = 10000
values = np.zeros(counts)   
for i in range (1, counts, 1):
    values[i] = random.random()

# the histogram of the data
n, bins, patches = plt.hist(values, 10, facecolor='green')

plt.xlabel('$x$')
plt.ylabel('Number of counts')
plt.title(r'Test of uniform distribution')
plt.axis([0, 1, 0, 1100])
plt.grid(True)

plt.show()

