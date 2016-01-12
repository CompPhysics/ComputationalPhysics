# -*- coding: utf-8 -*-

#This program demonstrates how to read command-line arguments,
#and how to read from- and write to file
#Translated to Python by Kyrre Ness Sjøbæk

import sys

if len(sys.argv) != 3:
    print "Usage: python " + sys.argv[0] + " infile outfile"
    sys.exit(1)
    
#We know there are three arguments - it should be safe to read them!
infile_name  = sys.argv[1];
outfile_name = sys.argv[2];

#Open the files
try:
    infile = open(infile_name,'r')
except:
    print "Oops! Could not read", infile_name
    sys.exit(1)

try:
    outfile = open(outfile_name, 'w')
except:
    print "Oops! Could not open", outfile_name, "for writing"
    sys.exit(1)

#Copy the contents of infile to outfile, line-by-line
for line in infile:
    outfile.write(line)

#Write a little son

#Close the files
infile.close()
outfile.close()
