# coding=utf-8
#Simple matrix inversion example
#Translated to Python by Kyrre Ness Sjøbæk

import numpy, computationalLib, math

myLib = computationalLib.pylib(inputcheck=False, cpp=False)

A = numpy.array([[1,3,4],[3,4,6],[4,6,8]],numpy.double)
print "Initial matrix A:\n", A

#LU-decompose A
matr = A.copy()
(matr,index,d) = myLib.luDecomp(matr);
print "LU-decomposed matrix:\n", matr
print "Index array:\n", index
print "Parity:\n", d

#Invert A column-by-column:
Ainv = numpy.zeros((3,3),numpy.double)
for j in xrange(3):
    col = numpy.zeros(3)
    col[j] = 1;

    Ainv[:,j] = myLib.luBackSubst(matr,index,col)

print "Inverse of A:\n", Ainv
    
#Test matrix inversion: AxAinv should be identity
test = True
ZERO = 1.0E-10
identity = numpy.zeros((3,3),numpy.double)
for i in xrange(3):
    for j in xrange(3):
        #Calculate and evaluate each element in resulting matrix:
        sum = 0.0;
        for k in xrange(3):
            sum = sum + A[i,k]*Ainv[k,j];
        if (i != j) and (sum > ZERO):
            test = False;
        if (i == j) and (math.fabs(sum - 1) > ZERO):
            test = False
        identity[i,j] = sum;
print "A*Ainv:\n", identity
if test:
    print "Status: OK"
else:
    print "Status: FAILED!"
            
