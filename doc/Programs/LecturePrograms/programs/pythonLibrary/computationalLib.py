# -*- coding: utf-8 -*-

import numpy, math
import os,sys


class pylib:
    """
    Implements many of the functions found in M.H.Jensens
    c++-library, used in Computational Physics. This is again heavily based
    on what is found in Numerical recipes.
    
    Ported to Python by
    Kyrre Ness Sjøbæk      (k DÅTT n DÅTT sjobak ÆTT fys DÅTT uio DÅTT no),
    Magnar Kopangen Bugge  (magnarkb ÆTT student DÅTT matnat DÅTT uio DÅTT no),
    Marit Sandstad         (marit DÅTT sandstad ÆTT fys DÅTT uio DÅTT no)
    """

    
    ZERO = 1.0E-10;
    """Used as a meassure of machine precision in some algos"""

    yTol = None
    """Vector or number giving the tolerance for the stepsize control functions"""
    guard = 0.95
    """Used for the adaptive stepsize control functions. Standard value (0.95) is just a guess..."""

    def __init__(self, inputcheck=False, cpp=False):
        """
        Constructor,
        Set inputcheck = True in order to do input checking (debug, slower)
        Set cpp = True to use compiled versions of the functions where aviable
        """
        self.inputcheck = inputcheck
        self.cpp        = cpp
        if cpp:
            sys.path.append(os.path.abspath(".") + "/cpp")
            global pylib_cpp
            import pylib_cpp

    def luDecomp(self, A):
        """
        LU-decomposes a matrix A, and returns the LU-decomposition of a
        rowwise permutation of A. Used in combination with luBackSubst function
        to solve an equation-set Ax=B.

        Returns: A tuple containing
        - LU-decomposed matrix (upper part and diagonal = matrix U, lower part = matrix L)
        - Array which records the row permutation from partial pivoting
        - Number which depends on the number of row interchanges was even (+1) or odd (-1)

        BIG FAT WARNING: Destroys input A in calling program!
        (A is set equal to the returned LU-decomposed matrix)
        Send it a copy if this is bad.

        This function has the ability to switch between Python and C++ backends, see __init__()
        """

        if self.inputcheck:
            self.checkSquare(A)

        d = 1.0; #Records row interchange (parity)
        N = A.shape[0]
        index = numpy.zeros(N,numpy.int32)

        if self.cpp:
            return self.luDecomp_cpp(A,N,index,d);
        else:
            return self.luDecomp_python(A,N,index,d);

    def luDecomp_cpp(self,A,N,index,d):
        """
        C++ backend for luDecomp, using routine
        in the library
        """
        pylib_cpp.ludcmp(A,index,d)
        return (A,index,d)

    def luDecomp_python(self,A,N,index,d):
        """
        Python backend for luDecomp
        """
        
        #Loop over rows to get scaling info, check if matrix is singular
        vv = numpy.zeros(N)
        for i in xrange(N):
            big = numpy.fabs(A[i]).max();
            if big < self.ZERO:
                raise SingularError(A)
            vv[i] = 1.0/big

        #Loop over columns, Crout's method
        for j in xrange(N):
            #i < j
            for i in xrange(j):
                sum = A.item(i,j)
                for k in xrange(i):
                    sum = sum - A.item(i,k)*A.item(k,j)
                A.itemset(i,j,sum);

            #i >= j
            imax = None;
            big = self.ZERO
            for i in xrange(j,N):
                sum = A.item(i,j)
                for k in xrange(j):
                    sum = sum - A.item(i,k)*A.item(k,j);
                A.itemset(i,j,sum);
                #Find biggest entry in vv (scaled by THIS sum)
                dum = vv[i]*math.fabs(sum)
                if dum >= big:
                    big = dum;
                    imax = i;
            #Do we need to interchange rows?
            if j != imax:
                #Oui!
                dum2 = A[imax].copy()
                A[imax] = A[j]
                A[j] = dum2
                d = d*(-1)
                vv[imax] = vv[j];
            index[j] = imax;

            #If the pivot element is to small, the matrix is singular
            #to working order. For some applications of singular matrices,
            #it is desirable to substitute in self.ZERO
            if math.fabs(A.item(j,j) < self.ZERO):
                A.itemset(j,j,self.ZERO)

            #Divide by pivot element
            if j < (N-1):
                A[(j+1):N,j] = A[(j+1):N,j]/A.item(j,j)

        return (A,index,d)
    
    def luBackSubst(self, A, index, b):
        """
        Back-substitution of LU-decomposed matrix
        Solves the set of linear equations A x = b of dimension n
        The input A is the LU-decomposed version of A obtained from
        pylib.luDecomp(),
        index is the pivoting permutation vector obtained
        from the same function, and
        b is the right-hand-side vector b as a numpy array.

        Returns the solution x as an numpy array.

        BIG FAT WARNING: Destroys input b in calling program!
        (b is set equal to x after calculation has finished)
        Send it a copy if this is bad.
        """

        if self.inputcheck:
            self.checkSquare(A)
            if A.shape[0] != b.shape[0]:
                raise SizeError("A:(%d,%d), b:(%d)" % (A.shape[0], A.shape[1], b.shape[0]))

        N = A.shape[0];
        
        if self.cpp:
            return self.luBackSubst_cpp(A,N,index,b)
        else:
            return self.luBackSubst_python(A,N,index,b)
    def luBackSubst_cpp(self, A, N, index, b):
        """
        C++ backend for luDecomp, using routine
        in the library
        """
        pylib_cpp.luBackSubst(A,index,b)
        return b
    
    def luBackSubst_python(self, A, N, index, b):
        """
        Python backend for luBackSubst
        """
        
        ii = -1
        for i in xrange(N):
            ip    = index[i]
            sum   = b[ip]
            b[ip] = b[i]
            if ii > -1:
                for j in xrange(ii,i):
                    sum = sum - A[i,j]*b[j]
            elif sum != 0:
                ii = i
            b[i] = sum;
        for i in xrange((N-1),-1,-1):
            sum = b[i]
            for j in xrange(i+1,N):
                sum = sum - A[i,j]*b[j]
            b[i] = sum/A[i,i]

        return b;
        
    def trapezoidal(self, a,b, n, func):
        """
        Integrate the function func
        using the trapezoidal rule from a to b,
        with n points. Returns value from numerical integration
        """
        step = (b-a)/float(n)
        
        sum = func(a)/float(2);
        for i in xrange(1,n):
            sum = sum + func(a+i*step)
        sum = sum + func(b)/float(2)
        
        return sum*step

    def simpson(self,a,b,n,func):
        """Same as trapezoidal, but use simpsons rule"""
        step = (b-a)/float(n)
        
        sum = func(a)/float(2);
        for i in xrange(1,n):
            sum = sum + func(a+i*step)*(3+(-1)**(i+1))
        sum = sum + func(b)/float(2)
        
        return sum*step/3.0

    def gausLegendre(self,a,b,N):
        """
        Used to calculate weights and meshpoints for
        gaus-legendre quadrature.

        Input:
        - a, b: Limits of integration
        - N:    Integration points

        Returns:
        - Numpy array of meshpoints
        - Numpy array of weights

        Method heavily inspired of gauleg() from Numerical recipes
        Note: This is not the same method as described in compendium;
        see NR for more info!
        """

        x = numpy.zeros(N)
        w = numpy.zeros(N)

        if self.cpp:
            return self.gausLegendre_cpp(a,b,x,w,N)
        else:
            return self.gausLegendre_python(a,b,x,w,N)

    def gausLegendre_python(self,a,b,x,w,N):
        """
        Python backend for gausLegendre
        """
        #Roots are symetric around the midpoint of the interval
        m  = (N+1)/2   #How many points above midpoint (including midpoint if N odd)
        xm = (b+a)/2.0 #Midpoint of integration interval
        xl = (b-a)/2.0 #Half-size of interval

        for i in xrange(m):
            #Find the i'th zero-point by Newtons method
            #Need an approximation to start from
            z = math.cos(math.pi*(i+0.75)/(N+0.5))
            inaccurate = True
            pp = 0.0
            while inaccurate:
                p1 = 1.0
                p2 = 0.0
                #Generate L_N(z) by recursion relation
                for j in xrange(N):
                    p3 = p2
                    p2 = p1
                    p1 = ((2.0*j+1.0)*z*p2-j*p3)/(j+1)
                #L_N(z) is now stored in p1
                #Generate derivative by use of standard relation
                pp = N*(z*p1-p2)/(z*z-1.0)
                z1=z
                z=z1-p1/pp
                if math.fabs(z-z1) < self.ZERO:
                    inaccurate = False
            #Scale the roots to the desired interval, and put in its symetric counterpart
            x[i]      = xm-xl*z;
            x[-(i+1)] = xm+xl*z
            #Compute the weight, symetric counterpart
            w[i]=2.0*xl/((1.0-z*z)*pp*pp)
            w[-(i+1)] = w[i] #small array, no need for temp to avoid seek
        return (x,w)

    def gausLegendre_cpp(self,a,b,x,w,N):
        """
        C++ backend for gausLegendre, using routine
        in the library
        """
        pylib_cpp.gausLegendre(a,b,x,w,N)
        return (x,w)
    
    def jacobi(self, A, inaccurate=False):
        """
        Computes all eigenvalues and eigenvectors of a real symetric matrix A,
        using Jacobi transformations. Note: This is the cyclic Jacobi method,
        not the original one. See Numerical recipes!

        Input:
         - A: Square symetric matrix which is to be diagonalized.
              Lower part (below diag) left untouched, upper part destroyed
        - inaccurate: optional bool value to decide the degree of accuracy
              of the method. default value is False. 

        Output:
         - d:    Vector containing the eigenvalues of A
         - v:    Matrix with columns representing the eigenvectors of A (normalized)
         - nrot: Number of Jacobi rotations required
              
        """

        if self.inputcheck:
            self.checkSquare(A)
            self.checkSymetric(A)
        
        n = len(A)
	v = numpy.eye(n)    #Initializing v to the identity matrix
	d = numpy.zeros(n)
        nrot = 0 	      # Counting the number of rotations
        if self.cpp:
            if inaccurate:
                print "Warning from pylib.jacobi: Inaccurate \
                ignored in cpp mode"
            return self.jacobi_cpp(A,d,v,n,nrot)
        
        else:
            return self.jacobi_python(A, d, v, n, nrot, inaccurate)
        
    def jacobi_cpp(self,A, d, v, n, nrot):
        """
        C++ backend for jacobi, using routine
        in the library
        """
        nrot = pylib_cpp.jacobi(A, d, v, n) #Wont return nrot through INOUT. Strange.
        return (d, v, nrot)
        
    def jacobi_python(self, A, d, v, n, nrot, inaccurate):
        """
        Python backend for jacobi
        """
	b = numpy.zeros(n)
	z = numpy.zeros(n)  #z will accumulate terms of the form t*a_pq

	for ip in xrange(n):

	    b[ip] = A[ip][ip] #Initializing b and d to be the
	    d[ip] = A[ip][ip] #diagonal elements of A

	
	for i in xrange(1,500):
	    sm = 0.0
	    for ip in xrange(n-1):
		for iq in xrange(ip + 1, n):
			sm += math.fabs(A[ip][iq]) #Summing magnitude of 
					      #off-diagonal elements
	    if sm == 0.0:
                return  d, v, nrot
            elif inaccurate and sm <= self.ZERO:
                return  d, v, nrot  
	    if i < 4:
		tresh = 0.2*sm/(n*n)
	    else:
		tresh = 0.0
	    for ip in xrange(n-1):
		for iq in xrange(ip + 1, n):
		    g = 100.0 * math.fabs(A[ip][iq])

                    if i > 4 and ((math.fabs(d[ip]) + g) == math.fabs(d[ip])) and ((math.fabs(d[iq]) + g) == math.fabs(d[iq])):
                        A[ip][iq] = 0.0
                    elif math.fabs(A[ip][iq]) > tresh:
                        h = d[iq] - d[ip]
                        if (math.fabs(h) + g) == math.fabs(h):
                            t = A[ip][iq]/h
                        else:
                            theta = 0.5 * h / A[ip][iq]
                            t = 1.0 / (math.fabs(theta) + math.sqrt(1.0 + theta*theta))
                            if theta < 0.0:
                                t = - t
                        c = 1.0/math.sqrt(1 + t * t)
                        s = t * c
                        tau = s/(1.0 + c)
                        h = t * A[ip][iq]
                        z[ip] -= h
                        z[iq] += h
                        d[ip] -= h
                        d[iq] += h
                        A[ip][iq] = 0.0
                        for j in xrange(0, ip):
                            self.jacobi_rot(A, s, tau, j, ip, j, iq)
                        for j in xrange(ip + 1, iq):
                            self.jacobi_rot(A, s, tau, ip, j, j, iq)                      
                        for j in xrange(iq + 1, n):
                            self.jacobi_rot(A, s, tau, ip, j, iq, j)
                        for j in xrange(0, n):
                            self.jacobi_rot(v, s, tau, j, ip, j, iq)
                        nrot += 1
            for ip in xrange(0, n):
                b[ip] += z[ip]
                d[ip] = b[ip]
                z[ip] = 0.0
        
	print "Nå blir du kastet ut! nrot = %d" %nrot
	sys.exit(0)

    def jacobi_rot(self, A, s, tau, i, j, k, l):
        """
        Sub-routine used by jacobi() to do one part of a rotation

        Input:
         - A:   Matrix to be sub-rotated
         - s:   sine of rotation angle phi
         - tau: Tau (helper quantity)
         - i,j,k,l: \"coordninates\" in matrix to rotate
        """
	g = A[i][j]
	h = A[k][l]
	A[i][j] = g - s * (h + g * tau)
	A[k][l] = h + s * (g - h * tau)
	
	return A

    def tred2(self, A):
        """
        Perform a Householder reduction of a real symetric matrix
        to tridiagonal form. See Numerical recipes for more info!

        Input:
        - A: NxN Numpy array containing the matrix to be tridiagonalized. Destroyed by tred2()!

        Output:
        - d: Numpy array length N containing the diagonal elements of the tridiagonalized matrix
        - e: Numpy array length N containing the off-diagonal elements of the tridiagonalized matrix

        - Input A is replaced by the orthogonal matrix effecting the transformation (NOT eigenvectors)
        => BIG FAT WARNING! A destroyed.
        """

        if self.inputcheck:
            self.checkSquare(A)
            self.checkSymetric(A)
            #A should be floating-point
            if type(A[0,0]) != numpy.float64:
                raise TypeError("Matrix not floating-point!")
            

        #Initialize variables
        N = A.shape[0]
        d = numpy.zeros(N)
        e = numpy.zeros(N)

        if self.cpp:
            return self.tred2_cpp(A,N,d,e)
        else:
            return self.tred2_python(A,N,d,e)

    def tred2_cpp(self,A,N,d,e):
        """
        C++ backend for tred2, using routine
        in the library
        """

        pylib_cpp.tred2(A,N,d,e);
        return (d,e)

    def tred2_python(self,A,N,d,e):
        """
        Python backend for tred2
        """
        #Loop backwards through the matrix (modified Housholder -- see NR)
        for i in xrange(N-1,0,-1): #Start, stop, step
            l = i-1 #Index of last entry in row before diagonal
            h=scale=0.0
            
            if (l > 0):
                #for k in xrange(l+1):
                #    scale += math.fabs(A.item(i,k))
                #Equivalent:
                scale = numpy.sum(numpy.fabs(A[i,:l+1]))
                if (scale == 0.0):
                    #Skip this transformation
                    e[i] = A[i,l]
                else:
                    for k in xrange(l+1):
                        #Scale the entries
                        A[i,k] /= scale
                        #Form sigma in h
                        h+=A[i,k]*A[i,k]
                    f = A[i,l]
                    if f >= 0.0:
                        g = - math.sqrt(h)
                    else:
                        g = math.sqrt(h)
                    e[i] = scale*g
                    h -= f*g
                    A[i,l] = f-g
                    f=0.0
                    for j in xrange(l+1):
                        A[j,i] = A[i,j]/h
                        g = 0.0
                        for k in xrange(k+1):
                            g += A[k,j]*A[i,k]
                        e[j] = g/h
                        f += e[j]*A[i,j]
                    hh = f/(h+h)
                    for j in xrange(l+1):
                        f = A[i,j]
                        e[j]=g = e[j] - hh*f;
                        for k in xrange(j+1):
                            A[j,k] -= (f*e[k] + g*A[i,k])

            else: #(l = 0 => i=1, last row (row 2 in math-notation)
                e[i]=A[i,l]

            d[i] = h;

        d[0] = e[0] = 0.0
        for i in xrange(N):
            l=i
            if (d[i] != 0.0): #Skipped when i=0
                for j in xrange (l):
                    g = 0.0
                    for k in xrange(l):
                        g+=A[i,k]*A[k,j]
                    for k in xrange(l):
                        A[k,j] -= g*A[k,i]
            d[i] = A[i,i]
            #Reset row/column to identity matrix for next iteration
            A[i,i] = 1.0
            for j in xrange(l):
                A[j,i]=A[i,j]=0.0

        return(d,e)
                    
                        
    def pythag(self,a,b):
        """
        Function which computes sqrt(a^2+b^2) without loss of precision. Used by tqli().
        """
        absa = abs(a)
        absb = abs(b)
        if absa > absb:
            return absa*math.sqrt(1.0+(absb/absa)**2)
        elif absb == 0:
            return 0.0
        else:
            return absb*math.sqrt(1.0+(absa/absb)**2)

    def sign(self,a,b):
        """
        Function which returns |a| * sgn(b).
        """
        if b < 0:
            return -abs(a)
        else:
            return abs(a)

    def tqli(self,d,e,z):
        """
        Function which finds the eigenvalues and eigenvectors of a tridiagonal symmetric
        matrix. This is a translation of the function tqli in lib.cpp. 

        Input:
        - d: diagonal elements of the matrix
        - e: off-diagonal elements of the matrix (first element is dummy)
        - z: unit matrix (if eigenvectors of an \"already tridiag\" matrix wanted),
             or matrix from tred2() if eigenvectors from a symetric matrix tridiagonalized by
             tred2() wanted.

        The elements of d after tqli() has finished are the eigenvalues of
        the tridiagonal matrix. The eigenvectors are stored in z.
        """

        if self.inputcheck:
            self.checkSquare(z)
            if len(d) != len(e) or len(d) != z.shape[0]:
                raise SizeError("len(d) = %d, len(e) = %d, z.shape[0] = %d" %(len(d),len(e),z.shape[0]))

        n = len(d) #n was given as an input parameter in the C++ library

        if self.cpp:
            self.tqli_cpp(d,e,n,z);
        else:
            self.tqli_python(d,e,n,z)
            
    def tqli_cpp(self,d,e,n,z):
        """
        C++ backend for tqli(), using routine
        in the library
        """
        pylib_cpp.tqli(d,e,n,z);

    def tqli_python(self,d,e,n,z):
        """
        Python backend for tqli()
        """
        for i in xrange(1,n):
            e[i-1] = e[i]
        e[n-1] = 0.0

        for l in xrange(n):
            iter = 0

            while True:

                m = l
                while m < n-1:
                    dd = abs(d[m]) + abs(d[m+1])
                    if abs(e[m]) + dd == dd:
                        break
                    m += 1
                    
                if m != l:
                    if iter == 30:
                        print '\nToo many iterations in tqli.\n'
                        return
                    iter += 1

                    g = (d[l+1] - d[l])/(2.0 * e[l])
                    r = self.pythag(g,1.0)
                    g = d[m]-d[l]+e[l]/(g+self.sign(r,g))
                    s = c = 1.0
                    p = 0.0
                                         
                    i = m - 1
                    while i >= l:
                        f = s * e[i]
                        b = c * e[i]
                        r = self.pythag(f,g)
                        e[i+1] = r
                        if r == 0.0:
                            d[i+1] -= p
                            e[m] = 0.0
                            break

                        s = f/r
                        c = g/r
                        g = d[i+1] - p
                        r = (d[i] - g) * s + 2.0 * c * b
                        p = s*r
                        d[i+1] = g + p
                        g = c * r - b

                        for k in xrange(n):
                            f = z[k,i+1]
                            z[k,i+1] = s * z[k,i] + c * f
                            z[k,i] = c * z[k,i] - s * f

                        i -= 1

                    if r == 0.0 and i >= l:
                        continue

                    d[l] -= p
                    e[l] = g
                    e[m] = 0.0
                else:
                    break


    def euler(self, y0, t0, te, N, deriv, filename=None):
        """
        General eulers method driver and stepper for
        N coupled differential eq's,
        fixed stepsize
        
        Input:
         - y0:       Vector containing initial values for y
         - t0:       Initial time
         - te:       Ending time
         - N:        Number of steps
         - deriv:    See rk4_step
         - filename: Optional, use if you want to write
                     data to file at each step.
                     Format used:
                     t y[0] y[1] ... (%10.15E)

        Output:
        If filename=None, return tuple containing:
         - time:  Array of times at which it has iterated over
         - yout:  N*len(y0) numpy array containing y for each timestep
        If filename specified, None is returned.
        """
        
        h = (te-t0)/float(N)
        t = t0;
        
        if filename == None:
            #Setup arrays
            time = numpy.zeros(N);
            yout = []
            #Inital values
            yout.append(y0);
            time[0] = t0;
            t = t0;
            
            #Loop over timesteps
            for i in xrange(1,N):
                yout.append(y + h*deriv(y,t));
                t = t0 + h*i;
                time[i] = t;
                
            return (time,yout)
        else:
            ofile = open(filename,'w')
            #Format string used for output file
            ostring = "%20.8E " + ("%20.8E "*len(y0)) + "\n"
            
            #Initial values
            y = y0
            t = t0
            
            foo = [t]; foo[1:] = y;
            ofile.write(ostring % tuple(foo))
        
            while (t < te):
                y = y + h*deriv(y,t)
                t +=h
                
                foo = [t]; foo[1:] = y;
                ofile.write(ostring % tuple(foo))
                
            ofile.close()
            return None

    def rk2(self, y0, t0, te, N, deriv, filename=None):
        """
        General RK2 driver for
        N coupled differential eq's,
        fixed stepsize
        
        Input:
         - y0:       Vector containing initial values for y
         - t0:       Initial time
         - te:       Ending time
         - N:        Number of steps
         - deriv:    See rk4_step
         - filename: Optional, use if you want to write
                     data to file at each step.
                     Format used:
                     t y[0] y[1] ... (%10.15E)

        Output:
        If filename=None, return tuple containing:
         - time:  Array of times at which it has iterated over
         - yout:  N*len(y0) numpy array containing y for each timestep
        If filename specified, None is returned.
        """
        
        h = (te-t0)/float(N)
        t = t0;
        
        if filename == None:
            #Setup arrays
            time = numpy.zeros(N);
            yout = []
            #Inital values
            yout.append(y0);
            time[0] = t0;
            t = t0;
            
            #Loop over timesteps
            for i in xrange(1,N):
                yout.append(self.rk2_step(yout[i-1],t,h,deriv));
                t = t0 + h*i;
                time[i] = t;
                
            return (time,yout)
        else:
            ofile = open(filename,'w')
            #Format string used for output file
            ostring = "%20.8E " + ("%20.8E "*len(y0)) + "\n"
            
            #Initial values
            y = y0
            t = t0
            
            foo = [t]; foo[1:] = y;
            ofile.write(ostring % tuple(foo))
        
            while (t < te):
                y = self.rk2_step(y,t,h,deriv)
                t +=h
                
                foo = [t]; foo[1:] = y;
                ofile.write(ostring % tuple(foo))
                
            ofile.close()
            return None
        
    def rk2_step(self,y,t,h,deriv):
        """
        General RK2 stepper for
        N coupled differential eq's

        Input:
            - y:      Array containing the y(t)
            - t:      Which time are we talking about?
            - h:      Stepsize
            - deriv:  Function that returns an array
                      containing dy/dt for each y at time t,
                      and takes as arguments an y-array, and time t.
        """

        k1 = h*deriv(y,t);
        k2 = h*deriv(y+k1/2.0,t+h/2.0)
        

        return y + k2;

    def rk4(self, y0, t0, te, N, deriv, filename=None):
        """
        General RK4 driver for
        N coupled differential eq's,
        fixed stepsize
        
        Input:
         - y0:       Vector containing initial values for y
         - t0:       Initial time
         - te:       Ending time
         - N:        Number of steps
         - deriv:    See rk4_step
         - filename: Optional, use if you want to write
                     data to file at each step.
                     Format used:
                     t y[0] y[1] ... (%10.15E)

        Output:
        If filename=None, return tuple containing:
         - time:  Array of times at which it has iterated over
         - yout:  N*len(y0) numpy array containing y for each timestep
        If filename specified, None is returned.
        """
        
        h = (te-t0)/float(N)
        t = t0;
        
        if filename == None:
            #Setup arrays
            time = numpy.zeros(N);
            yout = []
            #Inital values
            yout.append(y0);
            time[0] = t0;
            t = t0;
            
            #Loop over timesteps
            for i in xrange(1,N):
                yout.append(self.rk4_step(yout[i-1],t,h,deriv));
                t = t0 + h*i;
                time[i] = t;
                
            return (time,yout)
        else:
            ofile = open(filename,'w')
            #Format string used for output file
            ostring = "%20.8E " + ("%20.8E "*len(y0)) + "\n"
            
            #Initial values
            y = y0
            t = t0
            
            foo = [t]; foo[1:] = y;
            ofile.write(ostring % tuple(foo))
        
            while (t < te):
                y = self.rk4_step(y,t,h,deriv)
                t +=h
                
                foo = [t]; foo[1:] = y;
                ofile.write(ostring % tuple(foo))
                
            ofile.close()
            return None

    def rk4_step(self,y,t,h,deriv):
        """
        General RK4 stepper for
        N coupled differential eq's

        Input:
            - y:      Array containing the y(t)
            - t:      Which time are we talking about?
            - h:      Stepsize
            - deriv:  Function that returns an array
                containing dy/dt for each y at time t,
                and takes as arguments an y-array, and time t.
        """

        k1 = h*deriv(y,t);
        k2 = h*deriv(y+k1/2.0,t+h/2.0)
        k3 = h*deriv(y+k2/2.0,t+h/2.0)
        k4 = h*deriv(y+k3,t+h)

        return y + (k1 + 2*(k2+k3) + k4)/6.0

    def rk4Adaptive(self, y0, t0, te, h0, deriv, errorfunc, filename=None):
        """
        General RK4 driver for
        N coupled differential eq's,
        adaptive stepsize.

        Inspired by Numerical Recipes, and
        http://www.cofc.edu/lemesurierb/math545-2007/handouts/adaptive-runge-kutta.pdf
        plus some of my own extras.
        
        Input:
         - y0:        Vector containing initial values for y
         - t0:        Initial time
         - te:        Ending time
         - h0:        Initial guess for stepsize h
         - deriv:     See rk4_step()
         - errorfunc: Fuction that returns 0.0 if the step is accepted,
                      or else a returns the new stepsize
                      Expects input yOld, yNew, yErr, h, and t
         - filename:  Optional, use if you want to write
                      data to file at each step.
                      Format used:
                      t y[0] y[1] ... (%10.15E)

        Output:
        If filename=None, return tuple containing:
         - time:  Array of times at which it has iterated over
         - yout:  N*len(y0) numpy array containing y for each timestep
        If filename specified, None is returned.
        """

        # http://en.wikipedia.org/wiki/Pebkac
        if errorfunc == self.rk4Adaptive_stepsizeControl1 or errorfunc == self.rk4Adaptive_stepsizeControl2:
            if self.yTol == None:
                print "Please set yTol to use rk4Adaptive_stepsizeControl"
                sys.exit(0)

        if filename != None:
            ofile = open(filename,'w')
            #Format string used for output file
            ostring = "%20.8E " + ("%20.8E "*len(y0)) + "\n"
            foo = [t0]; foo[1:] = y0;
            ofile.write(ostring % tuple(foo))
        else:
            #Setup arrays
            time = []
            yout = []
            #Inital values
            yout.append(y0);
            time.append(t0);
            
        t = t0;
        i = 1
        h = h0
        yOld = y0
        
        #Loop until end of time
        while t < te:
            #Bigstep:
            y2h = self.rk4_step(yOld,t,2*h,deriv)
            #Doublestep
            yh  = self.rk4_step(yOld,t,h,deriv)
            yhh = self.rk4_step(yh,t+h,h,deriv)
            
            yErr = yhh - y2h;
            y    = yhh + yErr/15.0 # Estimate accurate to O(h^6)
            
            #Accept?
            newh = errorfunc(yOld,y,yErr,h,t)
            if not newh: #Accept!
                t +=2*h
                i += 1
                if filename != None:
                    foo = [t]; foo[1:] = y;
                    ofile.write(ostring % tuple(foo))
                else:
                    yout.append(y)
                    time.append(t)
                yOld = y
            else: #Change stepsize
                if newh > h: #To small step - accept and change
                    t +=2*h
                    i += 1
                    if filename != None:
                        foo = [t]; foo[1:] = y;
                        ofile.write(ostring % tuple(foo))
                    else:
                        yout.append(y)
                        time.append(t)
                    yOld = y
                    
                    h = newh
                else:
                    #To big step - reject and change
                    h = newh
        if filename != None:
            ofile.close()
            return None            
        else:
            return (time,yout)

    def rk4Adaptive_stepsizeControl1(self, yOld, yNew, yErr, h, t):
        """
        Standard stepsize control algo for adaptive RK4.
        This variety uses a fractional tolerance,
        usefull for most problems that don't cross zero.
        yTol should in this case be a number, and is interpreted
        as an accuracy requirement epsilon.

        Also see rk4Adaptive() and rk4Adaptive_stepsizeControl1()

        Please set the class variable pylib.yTol before use!
        """
        yTol = self.yTol*yNew
        hNew = self.guard*h*(numpy.fabs(yTol/yErr).min())**0.2
        if hNew < self.ZERO:
            print "**** h = ZERO at time %.7E - solution may be innacurate! ****" % lowh
            return self.ZERO
        else:
            return hNew
        
    def rk4Adaptive_stepsizeControl2(self, yOld, yNew, yErr, h, t):
        """
        Standard stepsize control algo for adaptive RK4.
        This variety uses a fixed tolerance, usefull for oscilatory
        problems. You should set yTol to something like
        epsilon*yMax, where yMax isthe largest values you get,
        and epsilon is the desired accuracy (say, 10^-6).

        Also see rk4Adaptive() and rk4Adaptive_stepsizeControl2()

        Please set the class variable pylib.yTol before use!
        """
        hNew = self.guard*h*(numpy.fabs(self.yTol/yErr).min())**0.2
        if hNew < self.ZERO:
            print "**** h = ZERO at time %.7E - solution may be innacurate! ****" % lowh
            return self.ZERO
        else:
            return hNew

    def checkSquare(self,A):
        """
        Checks if A is square, if not it raises exeption MatrixNotSquareError.
        Called by relevant methods if self.inputcheck=True
        """
        if A.shape[0] != A.shape[1]:
            raise MatrixNotSquareError(A)

    def checkSymetric(self,A):
        """
        Checks if A is symetric, if not it raises exception MatrixNotSymetricError.
        Called by relevant methods if self.inputcheck=True. Always call checkSquare first!
        """
        N = A.shape[0]
        delta = 0.0
        for i in xrange(1,N):
            for j in xrange(i,N):
                delta += abs(A[i,j]-A[j,i])
        if delta != 0.0:
            raise MatrixNotSymetricError(A,delta)

class MatrixError(Exception):
    """
    General class for (mathematical) errors in matrices passed.
    """
    def __init__(self,A):
        """
        Input:
        - A: The matrix which has an error
        """
        self.A = A;
    def __str__(self):
        return "Matrix = \n"+ str(self.A)

class MatrixNotSquareError(MatrixError):
    """
    Raised by a method requiring a square matrix if it is not square.
    """
    def __str__(self):
        return "Matrix not square, dimensions (%d,%d)" % (self.A.shape[0], self.A.shape[1])

class MatrixNotSymetricError(MatrixError):
    """
    Raised by a method requiring a symetric matrix if it is not symetric.
    Always raise MatrixNotSquareError before this one (a non-square matrix can't be symetric anyway).
    """
    def __init__(self,A,delta):
        """
        Input:
         - A: The matrix which has an error
         - delta: sum of abs(Aij-Aji)
        """        
        self.A = A;
        self.delta = delta
        
    def __str__(self):
        return "Matrix not symetric, delta = " + str(self.delta)

class SizeError(MatrixError):
    """
    Raised if a set of matrices/vectors that should have the same size
    have different sizes
    """
    def __init__(self, message):
        """
        Input:
         - message: String containing relevant info about the sizes and names
        """
        self.message = message
    def __str__(self):
        return "Matrices/vectors have different sizes: " + self.message

class SingularError(MatrixError):
    """
    Raised by numerical functions if a matrix is singular
    """
    def __str__(self):
        return "Matrix was singular:\n" + str(self.A)
