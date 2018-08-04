#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import random

# initialize the rng with a seed
random.seed()
# Hard coding of input parameters
Agents  = 100
MCcounts = 1000
Transactions = 10000
startMoney = 1.0
Lambda = 0.0
FinancialAgents = startMoney*np.ones(Agents)
for i in range (1, MCcounts, 1):
    for j in range (1, Transactions, 1):
        agent_i = int(Agents*random.random())
        agent_j = int(Agents*random.random())
        epsilon = random.random()
        if agent_i != agent_j:
           m1 = Lambda*FinancialAgents[agent_i] + (1-Lambda)*epsilon*    (FinancialAgents[agent_i] + FinancialAgents[agent_j])
           m2 = Lambda*FinancialAgents[agent_j] + (1-Lambda)*(1-epsilon)*(FinancialAgents[agent_i] + FinancialAgents[agent_j])
           FinancialAgents[agent_i] = m1
           FinancialAgents[agent_j] = m2

# the histogram of the data
n, bins, patches = plt.hist(FinancialAgents, 20, facecolor='green')

plt.xlabel('$x$')
plt.ylabel('Distribution of wealth')
plt.title(r'Money')
plt.axis([0, 10, 0, 100])
plt.grid(True)
plt.show()

