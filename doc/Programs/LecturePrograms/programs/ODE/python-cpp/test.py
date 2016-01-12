#!/usr/bin/env python

#This program shows same basic functionality as test.cpp

#Do stuff so we can import ./swig/pendelum.py
import os,sys
sys.path.append(os.path.abspath(".") + "/swig")
from pendelum import pendelum

testcase = pendelum()
testcase.initialise_keyboard();
testcase.euler();
testcase.euler_cromer();
testcase.midpoint();
testcase.euler_richardson();
testcase.half_step();
testcase.rk2();
testcase.rk4();
