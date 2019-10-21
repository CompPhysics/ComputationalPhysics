from  matplotlib import pyplot as plt
import numpy as np

# Define dimension of matrix and vectors
Dim = 4
#Setting up a transition probability matrix
TransitionMatrix = np.matrix('0.25 0.1111 0.375 0.3333; 0.5 0.2222 0.0 0.3333; 0.0 0.1111 0.375 0.0; 0.25 0.5556 0.25 0.3334')
# Making a copy of the transition matrix
W = TransitionMatrix
print(W)
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
for i in range(Dim):
    print(EigValues[i])
FifthEigvector = EigVectors[:,3]

#print(FifthEigvector)

difference = np.linalg.norm(FifthEigvector-wold,2)
print(difference)
eps = 1.0E-10
count = 0
while count < 10:
      for i in range(Dim):
          wnew[i] = np.dot(W[i,:],wold)
      count = count + 1
      print(count, wnew)
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
