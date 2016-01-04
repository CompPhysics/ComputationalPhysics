from  matplotlib import pyplot as plt
import numpy as np

# Define dimension of matrix and vectors
Dim = 2
#Setting up a transition probability matrix
TransitionMatrix = np.matrix('0.25 0.5; 0.75 0.5')
# Making a copy of the transition matrix
W = TransitionMatrix
print W
# our first state
wold = np.zeros(Dim)
wold[0] = 1.0
wnew = np.zeros(Dim)

# diagonalize and obtain eigenvalues, not necessarily sorted
EigValues, EigVectors = np.linalg.eig(TransitionMatrix)
# sort eigenvectors and eigenvalues
permute = EigValues.argsort()
EigValues = EigValues[permute]
EigVectors = EigVectors[:,permute]
for i in xrange(Dim):
    print EigValues[i]
FifthEigvector = EigVectors[:,1]

print FifthEigvector

difference = np.linalg.norm(FifthEigvector-wold,2)
print difference
eps = 1.0E-10
count = 0
while count < 10:
      for i in range(Dim):
          for j in range(Dim):
              wnew[i] = wnew[i]+W[i,j]*wold[j]
      count = count + 1
      print count, wnew
#      difference = np.linalg.norm(FifthEigvector-wnew,2)
#      print difference
      wold = wnew
#plt.plot(r, FirstEigvector**2 ,'b-',r, SecondEigvector**2 ,'g-',r, ThirdEigvector**2 ,'r-')
#plt.axis([0,4.6,0.0, 0.025])
#plt.xlabel(r'$r$')
#plt.ylabel(r'Radial probability $r^2|R(r)|^2$')
#plt.title(r'Radial probability distributions for three lowest-lying states')
#plt.savefig('eigenvector.pdf')
#plt.show()
