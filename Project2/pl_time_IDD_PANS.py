import numpy as np
import matplotlib.pyplot as plt
import sys
#from scipy.signal import welch, hanning
from scipy.signal import welch, hann
plt.rcParams.update({'font.size': 22})

# ***** read u
# /chalmers/users/lada/python-DES-code/channel-5200-IDD-PANS-inlet-synt-MTF271/u-time-history-i70.dat

plt.interactive(True)

n1=7500
n2=15000

data = np.genfromtxt("u-time-history-i70.dat", dtype=None)
t=data[n1:n2,0] #time
u1=data[n1:n2,1]   #v_1 at point 1
u2=data[n1:n2,2]   #v_1 at point 2
u3=data[n1:n2,3]   #v_1 at point 3
u4=data[n1:n2,4]   #v_1 at point 4
u5=data[n1:n2,5]   #v_1 at point 5

data = np.genfromtxt("w-time-history-i70.dat", dtype=None)
w1=data[n1:n2,1]   #v_3 at point 1
w2=data[n1:n2,2]   #v_3 at point 2
w3=data[n1:n2,3]   #v_3 at point 3
w4=data[n1:n2,4]   #v_3 at point 4
w5=data[n1:n2,5]   #v_3 at point 5

dx=0.1
dz=0.05

# plot time history 
fig1,ax1 = plt.subplots()
plt.subplots_adjust(left=0.20,bottom=0.20)

plt.plot(t,u1,'b-')
# plot every 30th time step
plt.plot(t[::30],u1[::30],'bo')
plt.plot(t,u3,'r-')
plt.xlabel('$t$')
plt.axis([10, 11, 9,25])
plt.ylabel('$u$')
plt.savefig('u-time.eps')

