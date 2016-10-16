# This python script reads from the command line the filename (its root) and
# the largest exponent of 10. This defines the number of mesh points.
# This script calls then an executable from a c++ or fortran code that solves
# a set of linear equations with a tridiagonal matrix defining the second derivative.
# It makes in turn the various plots as pdf files and finally sets up the basis for
# a report and its pertinent latex file

import sys, os
from  matplotlib import pyplot as plt
import numpy as np

# Command line arguments using sys.argv[]
try:
    filename = sys.argv[1]
    exponent = int(sys.argv[2])
except:
    print "Usage of this script", sys.argv[0], "infile", sys.argv[1], "Exponent", sys.argv[2]; sys.exit(1)

# Define command line text string 
cmdline  = './project1.x '+filename +' ' + str(exponent)
# Now run code, here c++ code  which has been compiled and linked
cmd = cmdline
failure = os.system(cmd)
if failure:
   print 'running project1 failed'; sys.exit(1)

# Start making figures looping over all exponents
for i in range(1,exponent+1):
#   define files to open data and make plots to
    fout = filename+str(i)
    figfile = fout+".pdf"
    data = np.loadtxt(fout)
    x = data[:,0]
    solution = data[:,1]
    exact = data[:,2]
    plt.axis([0,1,0, 1.0])
    numericalplot = plt.plot(x, solution, 'r:.', linewidth = 2.0, label = 'Numerical')
    exactplot = plt.plot(x, exact, 'm:v', linewidth = 2.0, label = 'Exact')
    plt.xlabel(r'$x$')
    plt.ylabel(r'Solutions')
    plt.savefig(figfile)
#   Then clean up
    plt.clf()


# Now prepare latex file, r in front avoids backslashes being treated
# as control chars in strings. What follows are plain  latex commands
preamb = r"""\documentclass[10pt,showpacs,preprintnumbers,footinbib,amsmath,amssymb,aps,prl,twocolumn,groupedaddress,superscriptaddress,showkeys]{revtex4-1}
\usepackage{graphicx}
\usepackage{dcolumn}
\usepackage{bm}
\usepackage[colorlinks=true,urlcolor=blue,citecolor=blue]{hyperref}
\usepackage{color}
\begin{document}
\title{Project 1}
\author{A.~N.~Author}
\affiliation{Department of Something, University of Somewhere, Outer Space}
\begin{abstract}
We present our Ferrari algorithm for solving linear equations. Our best algorithm runs as $4n$ FLOPS with $n$ the dimensionality of the matrix.
\end{abstract}
\maketitle

"""

figure = r"""\begin{figure}[hbtp]
\includegraphics[scale=0.4]{test1.pdf}
\caption{Exact and numerial solutions for $n=10$ mesh points.} 
\label{fig:n10points}
\end{figure}

"""


introduction = r"""\section{Introduction}

"""

theory = r"""\section{Theory, algorithms and methods}

"""

results = r"""\section{Results and discussions}

"""

conclusions = r"""\section{Conclusions}

"""

references = r"""\begin{thebibliography}{99}
\bibitem{miller2006} G.~A.~Miller, A.~K.~Opper, and E.~J.~Stephenson, Annu.~Rev.~Nucl.~Sci.~{\bf 56}, 253 (2006).
\end{thebibliography}

"""

# Dump to file:
filename = 'ReportProject1'
f = file(filename + '.tex', "w")
f.write(preamb)
f.write(introduction)
f.write(theory)
f.write(results)
f.write(figure)
f.write(conclusions)
f.write(references)
f.write("""\end{document}""")
f.close()









