import numpy as np
from math import *
import matplotlib.pyplot as plt

def solver(m, dt, t0):
    """Solve the difference equations for H and L over m years
    with time step dt (measured in years."""

    num_intervals = int(m/float(dt))
    print num_intervals
    t = np.linspace(t0, t0 + m, num_intervals+1)
    x = np.zeros(t.size)
    y = np.zeros(t.size)
    vx = np.zeros(t.size)
    vy = np.zeros(t.size)
    r = np.zeros(t.size)
    v = np.zeros(t.size)
    x[0] = 1.0
    y[0] = 0.0
    vx[0] = 2.0*pi
    vy[0] = 0.0
    pi4 = 4.0*pi*pi
    for n in range(0, len(t)-1):
        x[n+1] = x[n] + dt*vx[n]
        y[n+1] = y[n] + dt*vy[n]
        r3 = (x[n]*x[n]+y[n]*y[n])**3
        vx[n+1] = vx[n] -dt*pi4*x[n]/r3
        vy[n+1] = vy[n] -dt*pi4*y[n]/r3
        v[n+1] = sqrt(vx[n+1]*vx[n+1]+vy[n+1]*vy[n+1])
        r[n+1] = sqrt(x[n+1]*x[n+1]+y[n+1]*y[n+1])

    return r, v, t
# Simulate using the model
m =20 # 20 years
dt =0.01 # stepsize
t0 =0.0
r, v, t = solver(m, dt,t0)
# Visualize simulations and data
plt.plot(t, r, 'b-+')
plt.xlabel('r')
plt.ylabel('velocity')
plt.axis([0, 20.0, 0.0, 1.0])
plt.title(r'Velocity versus position')
plt.savefig('SunEarth.pdf')
plt.show()
