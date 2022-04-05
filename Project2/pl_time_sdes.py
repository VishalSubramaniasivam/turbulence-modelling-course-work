import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch, hann
from IPython import display
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
#plt.rcParams.update({'font.size': 20})

# makes sure figures are updated when using ipython
display.clear_output()

# ***** read u
data = np.genfromtxt("u_w_time_5nodes_sdes.dat", dtype=None)

u1=data[:,0]   #v_1 at point 1
u2=data[:,1]   #v_1 at point 2
u3=data[:,2]   #v_1 at point 3
u4=data[:,3]   #v_1 at point 4
u5=data[:,4]   #v_1 at point 5

w1=data[:,5]   #w_1 at point 1
w2=data[:,6]   #w_1 at point 2
w3=data[:,7]   #w_1 at point 3
w4=data[:,8]   #w_1 at point 4
w5=data[:,9]   #w_1 at point 5

#print("u1=",u1)

dx=3.2/32
dt= 0.25*dx/20
t_tot=dt*len(u1)

t = np.linspace(0,t_tot,len(u1))

# %%%%%%%%%%%%%%%% plotting section %%%%%%%%%%%%%%%%%%%%%%%%%%
# plot u
#fig1 = plt.figure("Figure 1")
#plt.plot(t,u1,'b--')
#plt.plot(t,u4,'r-')
#plt.xlabel("t")
#plt.ylabel("u")
#plt.savefig('utime_python.eps',bbox_inches='tight')

####################### # zoom
#fig2 = plt.figure("Figure 2")
#plt.plot(t,u1,'b--')
#plt.plot(t,u4,'r-')
#plt.xlabel("t")
#plt.ylabel("u")
#plt.axis([6, 7, 10, 22])
#plt.savefig('utime_zoom_python.eps',bbox_inches='tight')

# %% Assignment S1 - Time History
# Plot Time Variation of V1 at all nodes
fig1,axes = plt.subplots(2,1,constrained_layout=True)
axes = axes.reshape(-1)
axes[0].plot(t,u1,color='red',label=r'$y/\delta$ = 0.0028 (Node 1)')
axes[0].plot(t,u2,color='orange',label=r'$y/\delta$ = 0.0203 (Node 2)')
axes[0].plot(t,u3,color='blue',label=r'$y/\delta$ = 0.0364 (Node 3)')
axes[0].plot(t,u4,color='green',label=r'$y/\delta$ = 0.0645 (Node 4)')
axes[0].plot(t,u5,color='black',label=r'$y/\delta$ = 0.2 (Node 5)')
axes[1].plot(t,u1,color='red',label=r'$y/\delta$ = 0.0028 (Node 1)')
axes[1].plot(t,u2,color='orange',label=r'$y/\delta$ = 0.0203 (Node 2)')
axes[1].plot(t,u3,color='blue',label=r'$y/\delta$ = 0.0364 (Node 3)')
axes[1].plot(t,u4,color='green',label=r'$y/\delta$ = 0.0645 (Node 4)')
axes[1].plot(t,u5,color='black',label=r'$y/\delta$ = 0.2 (Node 5)')
axes[1].set_xlim(6,7)
axes[1].set_ylim(10,27)
for ax in axes:
    ax.set_xlabel(r'Time ($t$)')
    ax.set_ylabel(r'$\overline{u}/u_{\tau}$')
    ax.legend(fontsize=16,loc='best')
    ax.tick_params(labelsize='medium')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True, linestyle='-.')
fig1.suptitle(r'Time variation of $\overline{u}$ at various nodes',fontsize=22)

# Plot Time Variation of V3 at all nodes
fig2,axes = plt.subplots(2,1,constrained_layout=True)
axes = axes.reshape(-1)
axes[0].plot(t,w1,color='red',label=r'$y/\delta$ = 0.0028 (Node 1)')
axes[0].plot(t,w2,color='orange',label=r'$y/\delta$ = 0.0203 (Node 2)')
axes[0].plot(t,w3,color='blue',label=r'$y/\delta$ = 0.0364 (Node 3)')
axes[0].plot(t,w4,color='green',label=r'$y/\delta$ = 0.0645 (Node 4)')
axes[0].plot(t,w5,color='black',label=r'$y/\delta$ = 0.2 (Node 5)')
axes[1].plot(t,w1,color='red',label=r'$y/\delta$ = 0.0028 (Node 1)')
axes[1].plot(t,w2,color='orange',label=r'$y/\delta$ = 0.0203 (Node 2)')
axes[1].plot(t,w3,color='blue',label=r'$y/\delta$ = 0.0364 (Node 3)')
axes[1].plot(t,w4,color='green',label=r'$y/\delta$ = 0.0645 (Node 4)')
axes[1].plot(t,w5,color='black',label=r'$y/\delta$ = 0.2 (Node 5)')
axes[1].set_xlim(6,7)
axes[1].set_ylim(-3,2)
for ax in axes:
    ax.set_xlabel(r'Time ($t$)')
    ax.set_ylabel(r'$\overline{w}/u_{\tau}$')
    ax.legend(fontsize=16,loc='best')
    ax.tick_params(labelsize='medium')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True, linestyle='-.')
fig2.suptitle(r'Time variation of $\overline{w}$ at various nodes',fontsize=22)

# Compute Auto-Corelation
u1_fluct = u1 - np.mean(u1)
u2_fluct = u2 - np.mean(u2)
u3_fluct = u3 - np.mean(u3)
u4_fluct = u4 - np.mean(u4)
u5_fluct = u5 - np.mean(u5)

B11_1 = np.correlate(u1_fluct,u1_fluct,'full')
nmax = np.argmax(B11_1)
magmax = np.max(B11_1)
norm_B11_1 = B11_1[nmax:]/magmax

B11_2 = np.correlate(u2_fluct,u2_fluct,'full')
nmax = np.argmax(B11_2)
magmax = np.max(B11_2)
norm_B11_2 = B11_2[nmax:]/magmax

B11_3 = np.correlate(u3_fluct,u3_fluct,'full')
nmax = np.argmax(B11_3)
magmax = np.max(B11_3)
norm_B11_3 = B11_3[nmax:]/magmax

B11_4 = np.correlate(u4_fluct,u4_fluct,'full')
nmax = np.argmax(B11_4)
magmax = np.max(B11_4)
norm_B11_4 = B11_4[nmax:]/magmax

B11_5 = np.correlate(u5_fluct,u5_fluct,'full')
nmax = np.argmax(B11_5)
magmax = np.max(B11_5)
norm_B11_5 = B11_5[nmax:]/magmax

# Plotting Auto-Corelation
imax = 2000
fig3,axes = plt.subplots(2,1,constrained_layout=True)
axes[0].plot(t[0:imax],norm_B11_1[0:imax],color='red',label=r'$y/\delta$ = 0.0028 (Node 1)')
axes[0].plot(t[0:imax],norm_B11_2[0:imax],color='orange',label=r'$y/\delta$ = 0.0203 (Node 2)')
axes[0].plot(t[0:imax],norm_B11_3[0:imax],color='blue',label=r'$y/\delta$ = 0.0364 (Node 3)')
axes[1].plot(t[0:imax],norm_B11_4[0:imax],color='cyan',label=r'$y/\delta$ = 0.0645 (Node 4)')
axes[1].plot(t[0:imax],norm_B11_5[0:imax],color='black',label=r'$y/\delta$ = 0.2 (Node 5)')
for ax in axes:
    ax.set_xlabel(r'Time ($t$)')
    ax.set_ylabel(r'$B_{uu}$')
    ax.legend(fontsize=16,loc='best')
    ax.tick_params(labelsize='medium')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True, linestyle='-.')
    ax.set_xlim([0,0.07])
fig3.suptitle(r'Auto-Corelation ($B_{uu}$) at various nodes',fontsize=22)

# Compute Integral Timescale
# All normalised auto-corelations go negative. Find the instance when it goes
# negative and integrate till that point..
norm_B11_1_M = norm_B11_1[0:np.argmax(norm_B11_1<=0)]
norm_B11_2_M = norm_B11_2[0:np.argmax(norm_B11_2<=0)]
norm_B11_3_M = norm_B11_3[0:np.argmax(norm_B11_3<=0)]
norm_B11_4_M = norm_B11_4[0:np.argmax(norm_B11_4<=0)]
norm_B11_5_M = norm_B11_5[0:62]
#norm_B11_5_M = norm_B11_5[0:np.argmax(norm_B11_5<=0)]

int_T_1 = np.trapz(norm_B11_1_M) * dt
int_T_2 = np.trapz(norm_B11_2_M) * dt
int_T_3 = np.trapz(norm_B11_3_M) * dt
int_T_4 = np.trapz(norm_B11_4_M) * dt
int_T_5 = np.trapz(norm_B11_5_M) * dt