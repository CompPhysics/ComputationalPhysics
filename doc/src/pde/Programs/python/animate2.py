import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time, glob, os

def f(x, m, s):
    return (1.0/(np.sqrt(2*np.pi)*s))*np.exp(-0.5*((x-m)/s)**2)

m = 0
s_max = 2
s_min = 0.2
x = np.linspace(m -3*s_max, m + 3*s_max, 1000)
s_values = np.linspace(s_max, s_min, 30)
# f is max for x=m; smaller s gives larger max value
max_f = f(m, m, s_min)

# Make a first plot (save the lines objects returned from plt.plot)
fig = plt.figure()
plt.axis([x[0], x[-1], -0.1, max_f])
lines = plt.plot([], [])
plt.xlabel('x')
plt.ylabel('f')

# Function to return the background plot in the animation
def init():
    lines[0].set_data([], [])  # empty plot
    return lines

# Function to return a frame in the movie
def frame(args):
    frame_no, s, x, lines = args
    y = f(x, m, s)
    lines[0].set_data(x, y)
    # Does not work: lines[0].set_label('s=%4.2f' % s)
    # Does not work: plt.legend(['s=%4.2f' % s])
    # Does not work: plt.savefig('tmp_%04d.png' % frame_no)
    return lines

# Construct list of all arguments to frame function
# (each call sends frame number, s value, x array, and lines list)
all_args = [(frame_no, s, x, lines)
            for frame_no, s in enumerate(s_values)]

# Run the animation
anim = animation.FuncAnimation(
    fig, frame, all_args, interval=150, init_func=init, blit=True)

# Make movie file in MP4 format
anim.save('movie1.mp4', fps=5)
plt.show()
