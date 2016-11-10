import numpy as np
from  matplotlib import pyplot as plt
# Load in data file
data = np.loadtxt("autocor.dat")
data1 = np.loadtxt("automersenne.dat")
# Make arrays containing x-axis and binding energies as function of A
x = data[:,0]
corr = data[:,1]
corr2 = data1[:,1]
plt.plot(x, corr ,'ro', x, corr2, 'b')
plt.axis([0,1000,-0.2, 1.1])
plt.xlabel(r'$d$')
plt.ylabel(r'$C_d$')
plt.title(r'autocorrelation function for RNG')
plt.savefig('autocorr.pdf')
plt.show()
