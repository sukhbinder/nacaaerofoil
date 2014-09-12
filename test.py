"""
Testing aefoilnaca

"""

import aefoilnaca
import matplotlib.pyplot as plt

c = 10.0 # Chord length
t=0.15   # the maximum relative thickness
n=80     # number of point
m=0.02   # the maximum camber
p=0.4    # the location of maximum camber


y = aefoilnaca.naca4_sym(c,t,n)
yy = aefoilnaca.naca4_cam(0.02,0.4,t,c,80)

plt.plot(y[0,:],y[1,:])
plt.plot(yy[0,:],yy[1,:])
plt.show()
