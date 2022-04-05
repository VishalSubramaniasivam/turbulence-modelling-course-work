import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from dphidx_dy import dphidx_dy
from IPython import display
import numpy.ma as ma
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
plt.rcParams.update({'font.size': 22})

# makes sure figures are updated when using ipython
display.clear_output()


re =9.36e+5
viscos =1/re

xy_hump_fine = np.loadtxt("xy_hump_fine.dat")
x1=xy_hump_fine[:,0]
y1=xy_hump_fine[:,1]


nim1=int(x1[0])
njm1=int(y1[0])

ni=nim1+1
nj=njm1+1


x=x1[1:]
y=y1[1:]

x_2d=np.reshape(x,(njm1,nim1))
y_2d=np.reshape(y,(njm1,nim1))

x_2d=np.transpose(x_2d)
y_2d=np.transpose(y_2d)

# compute cell centers
xp2d= np.zeros((ni,nj))
yp2d= np.zeros((ni,nj))

for jj in range (0,nj):
   for ii in range (0,ni):

      im1=max(ii-1,0)
      jm1=max(jj-1,0)

      i=min(ii,nim1-1)
      j=min(jj,njm1-1)


      xp2d[ii,jj]=0.25*(x_2d[i,j]+x_2d[im1,j]+x_2d[i,jm1]+x_2d[im1,jm1])
      yp2d[ii,jj]=0.25*(y_2d[i,j]+y_2d[im1,j]+y_2d[i,jm1]+y_2d[im1,jm1])

# read data file
vectz = np.loadtxt("vectz_fine.dat")
ntstep=vectz[0]
ni=int(vectz[1])
nj=int(vectz[2])
nk=int(vectz[3])
n=len(vectz)

#            write(48,*)uvec(i,j)
#            write(48,*)vvec(i,j)
#            write(48,*)fk2d(i,j)
#            write(48,*)uvec2(i,j)
#            write(48,*)vvec2(i,j)
#            write(48,*)wvec2(i,j)
#            write(48,*)uvvec(i,j)
#            write(48,*)p2d(i,j)
#            write(48,*)rk2d(i,j)
#            write(48,*)vis2d(i,j)
#            write(48,*)dissp2d(i,j)
#            write(48,*)wvec(i,j)
#            write(48,*)vtvec(i,j)
#            write(48,*)tvec(i,j)


nn=14
nst=3
iu=range(nst+1,n,nn)
iv=range(nst+2,n,nn)
ifk=range(nst+3,n,nn)
iuu=range(nst+4,n,nn)
ivv=range(nst+5,n,nn)
iww=range(nst+6,n,nn)
iuv=range(nst+7,n,nn)
ip=range(nst+8,n,nn)
ik=range(nst+9,n,nn)
ivis=range(nst+10,n,nn)
idiss=range(nst+11,n,nn)

u=vectz[iu]/ntstep
v=vectz[iv]/ntstep
fk=vectz[ifk]/ntstep
uu=vectz[iuu]/ntstep
vv=vectz[ivv]/ntstep
ww=vectz[iww]/ntstep
uv=vectz[iuv]/ntstep
p=vectz[ip]/ntstep
k_model=vectz[ik]/ntstep
vis=vectz[ivis]/ntstep
diss=vectz[idiss]/ntstep

# uu is total inst. velocity squared. Hence the resolved turbulent resolved stresses are obtained as
uu=uu-u**2
vv=vv-v**2
uv=uv-u*v

p_2d=np.reshape(p,(ni,nj))
u_2d=np.reshape(u,(ni,nj))
v_2d=np.reshape(v,(ni,nj))
fk_2d=np.reshape(fk,(ni,nj))
uu_2d=np.reshape(uu,(ni,nj))
uv_2d=np.reshape(uv,(ni,nj))
vv_2d=np.reshape(vv,(ni,nj))
ww_2d=np.reshape(ww,(ni,nj))
k_model_2d=np.reshape(k_model,(ni,nj))
vis_2d=np.reshape(vis,(ni,nj)) #this is to total viscosity, i.e. vis_tot=vis+vis_turb
diss_2d=np.reshape(diss,(ni,nj)) #this is to total viscosity, i.e. vis_tot=vis+vis_turb

# set fk_2d=1 at upper boundary
fk_2d[:,nj-1]=fk_2d[:,nj-2]

# Plot Mesh for Visualization
a = np.transpose(xp2d)
b = np.transpose(yp2d)
c = np.transpose(x_2d)
d = np.transpose(y_2d)
fig1a = plt.figure()
axes = fig1a.add_axes([0.1,0.1,0.8,0.8])
#axes.plot(xp2d,yp2d,'r.')
#axes.plot(a,b,'r.')
axes.plot(x_2d,y_2d,'k-')
axes.plot(c,d,'k-')
#axes.plot(x_2d,y_2d,'b.')
#axes.plot(c,d,'b.')
axes.tick_params(axis="x", labelsize=15)
axes.tick_params(axis="y", labelsize=15)
axes.set_xlabel('x [m]',fontsize=20)
axes.set_ylabel('y [m]',fontsize=20)
axes.set_title('Computational mesh in the domain',fontsize=24)

dz=0.2/32
dx = np.zeros([ni,nj])
dy = np.zeros([ni,nj])
dx_dummy = x_2d[1:,1:] - x_2d[:-1,1:]
dy_dummy = y_2d[1:,1:] - y_2d[1:,:-1]
dx[1:-1,1:-1] = dx_dummy
dy[1:-1,1:-1] = dy_dummy
dx[0,:] = dx[1,:]
dx[:,0] = dx[:,1]
dx[-1,:] = dx[-2,:]
dx[:,-1] = dx[:,-2]
dy[0,:] = dy[1,:]
dy[:,0] = dy[:,1]
dy[-1,:] = dy[-2,:]
dy[:,-1] = dy[:,-2]
delta = (dx*dy*dx)**(1/3)
#delta = np.maximum.reduce([dx,dy,dz*np.ones([ni,nj])])

x065_off=np.genfromtxt("x065_off.dat",dtype=None,comments="%")
x066_off=np.genfromtxt("x066_off.dat",dtype=None,comments="%")
x080_off=np.genfromtxt("x080_off.dat",dtype=None,comments="%")
x090_off=np.genfromtxt("x090_off.dat",dtype=None,comments="%")
x100_off=np.genfromtxt("x100_off.dat",dtype=None,comments="%")
x110_off=np.genfromtxt("x110_off.dat",dtype=None,comments="%")
x120_off=np.genfromtxt("x120_off.dat",dtype=None,comments="%")
x130_off=np.genfromtxt("x130_off.dat",dtype=None,comments="%")


# compute the gradient
dudx,dudy=dphidx_dy(x_2d,y_2d,u_2d)


#*************************
# plot u
fig1b,ax1 = plt.subplots()
xx=0.65;
i1 = (np.abs(xx-xp2d[:,1])).argmin()  # find index which closest fits xx
plt.plot(u_2d[i1,:],yp2d[i1,:],'b-')
plt.plot(x065_off[:,2],x065_off[:,1],'bo')
plt.xlabel("$U$")
plt.ylabel("$y$")
plt.title("$x=0.65$")
plt.axis([0, 1.3,0.115,0.2])
# Create inset of width 30% and height 40% of the parent axes' bounding box
# at the lower left corner (loc=3)
# upper left corner (loc=2)
# use borderpad=1, i.e.
# 22 points padding (as 22pt is the default fontsize) to the parent axes
axins1 = inset_axes(ax1, width="40%", height="30%", loc=2, borderpad=1)
plt.plot(u_2d[i1,:],yp2d[i1,:],'b-')
plt.axis([0, 1.3,0.115,0.13])
# reduce fotnsize 
axins1.tick_params(axis = 'both', which = 'major', labelsize = 10)
# Turn ticklabels of insets off
axins1.tick_params(labelleft=False, labelbottom=False)
plt.plot(x065_off[:,2],x065_off[:,1],'bo')
#plt.savefig('u065_hump_python.eps',bbox_inches='tight')

# %% Assignment T4 - Location of Interface
# RANS to LES where fk goes down from 1 to 0.4
C_mu = 0.09
k_res_2d = 0.5*(uu_2d+vv_2d+ww_2d)
k_tot_2d = k_model_2d + k_res_2d
fk_def_2d = k_model_2d / k_tot_2d
Lt = (k_tot_2d**1.5)/diss_2d
fk_ref_2d = (C_mu**-0.5) * (delta/Lt)**(2/3)

switch_i = np.argmax((fk_2d[1:-1,1:-1]<=0.4),axis=1) + 1
switch_wd = np.zeros(len(switch_i))
for i in range (nim1-1):
    switch_wd[i] = yp2d[i+1,switch_i[i]] - yp2d[i+1,0]

fig2,axes = plt.subplots(1,1,constrained_layout=True)
axes.plot(xp2d[1:-1,0],switch_wd,label='D-PANS model')
axes.legend(fontsize=16,loc='best')
axes.xaxis.set_minor_locator(AutoMinorLocator())
axes.yaxis.set_minor_locator(AutoMinorLocator())
axes.tick_params(labelsize='medium')
axes.grid(True, linestyle='-.')
axes.set_xlabel('x (m)')
axes.set_ylabel(r'$y-y_{wall}$ (m)')
fig2.suptitle(r"Location of the switch (from RANS to LES) wrt the lower wall",fontsize=22)

xx = np.array([0,0.65,0.8,0.9,1,1.1,1.2,1.3,2])
i1 = list(map(lambda m : (np.abs(m-xp2d[:,1])).argmin(),xx))

fig3,axes = plt.subplots(2,4,constrained_layout=True)
axes = axes.reshape(-1)
axes[0].plot(fk_2d[i1[0],1:],yp2d[i1[0],1:]-yp2d[i1[0],0],'b-o',ms=4,label=r'$f_{k}$ in D-PANS')
axes[1].plot(fk_2d[i1[1],1:],yp2d[i1[1],1:]-yp2d[i1[1],0],'b-o',ms=4,label=r'$f_{k}$ in D-PANS')
axes[2].plot(fk_2d[i1[2],1:],yp2d[i1[2],1:]-yp2d[i1[2],0],'b-o',ms=4,label=r'$f_{k}$ in D-PANS')
axes[3].plot(fk_2d[i1[3],1:],yp2d[i1[3],1:]-yp2d[i1[3],0],'b-o',ms=4,label=r'$f_{k}$ in D-PANS')
axes[4].plot(fk_2d[i1[4],1:],yp2d[i1[4],1:]-yp2d[i1[4],0],'b-o',ms=4,label=r'$f_{k}$ in D-PANS')
axes[5].plot(fk_2d[i1[5],1:],yp2d[i1[5],1:]-yp2d[i1[5],0],'b-o',ms=4,label=r'$f_{k}$ in D-PANS')
axes[6].plot(fk_2d[i1[6],1:],yp2d[i1[6],1:]-yp2d[i1[6],0],'b-o',ms=4,label=r'$f_{k}$ in D-PANS')
axes[7].plot(fk_2d[i1[7],1:],yp2d[i1[7],1:]-yp2d[i1[7],0],'b-o',ms=4,label=r'$f_{k}$ in D-PANS')
axes[0].plot(fk_def_2d[i1[0],1:],yp2d[i1[0],1:]-yp2d[i1[0],0],'r-o',ms=4,label=r'$f_{k} = k/k_{tot}$')
axes[1].plot(fk_def_2d[i1[1],1:],yp2d[i1[1],1:]-yp2d[i1[1],0],'r-o',ms=4,label=r'$f_{k} = k/k_{tot}$')
axes[2].plot(fk_def_2d[i1[2],1:],yp2d[i1[2],1:]-yp2d[i1[2],0],'r-o',ms=4,label=r'$f_{k} = k/k_{tot}$')
axes[3].plot(fk_def_2d[i1[3],1:],yp2d[i1[3],1:]-yp2d[i1[3],0],'r-o',ms=4,label=r'$f_{k} = k/k_{tot}$')
axes[4].plot(fk_def_2d[i1[4],1:],yp2d[i1[4],1:]-yp2d[i1[4],0],'r-o',ms=4,label=r'$f_{k} = k/k_{tot}$')
axes[5].plot(fk_def_2d[i1[5],1:],yp2d[i1[5],1:]-yp2d[i1[5],0],'r-o',ms=4,label=r'$f_{k} = k/k_{tot}$')
axes[6].plot(fk_def_2d[i1[6],1:],yp2d[i1[6],1:]-yp2d[i1[6],0],'r-o',ms=4,label=r'$f_{k} = k/k_{tot}$')
axes[7].plot(fk_def_2d[i1[7],1:],yp2d[i1[7],1:]-yp2d[i1[7],0],'r-o',ms=4,label=r'$f_{k} = k/k_{tot}$')
for ax,xloc in zip(axes,xx):
    ax.plot(0.4*np.ones(20),np.linspace(0,0.2,20)-yp2d[i1[0],0],'k--',label=r'$f_{k}$ = 0.4')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(title='At X = '+str(xloc)+' m',title_fontsize=14,fontsize=14,loc='best')
    ax.tick_params(labelsize='medium')
    ax.grid(True, linestyle='-.')
    ax.set_ylim([0,0.2])
    if (ax==axes[0] or ax==axes[4]):
        ax.set_ylabel("$y-y_{wall}$")
    if (ax==axes[4] or ax==axes[5] or ax==axes[6] or ax==axes[7]):
        ax.set_xlabel("$f_{k}$")
fig3.suptitle(r"$f_{k}$ computed at various locations",fontsize=22)

fig4,axes = plt.subplots(2,4,constrained_layout=True)
axes = axes.reshape(-1)
axes[0].plot(fk_2d[i1[0],1:],yp2d[i1[0],1:]-yp2d[i1[0],0],'b-o',ms=4,label=r'$f_{k}$ in D-PANS')
axes[1].plot(fk_2d[i1[1],1:],yp2d[i1[1],1:]-yp2d[i1[1],0],'b-o',ms=4,label=r'$f_{k}$ in D-PANS')
axes[2].plot(fk_2d[i1[2],1:],yp2d[i1[2],1:]-yp2d[i1[2],0],'b-o',ms=4,label=r'$f_{k}$ in D-PANS')
axes[3].plot(fk_2d[i1[3],1:],yp2d[i1[3],1:]-yp2d[i1[3],0],'b-o',ms=4,label=r'$f_{k}$ in D-PANS')
axes[4].plot(fk_2d[i1[4],1:],yp2d[i1[4],1:]-yp2d[i1[4],0],'b-o',ms=4,label=r'$f_{k}$ in D-PANS')
axes[5].plot(fk_2d[i1[5],1:],yp2d[i1[5],1:]-yp2d[i1[5],0],'b-o',ms=4,label=r'$f_{k}$ in D-PANS')
axes[6].plot(fk_2d[i1[6],1:],yp2d[i1[6],1:]-yp2d[i1[6],0],'b-o',ms=4,label=r'$f_{k}$ in D-PANS')
axes[7].plot(fk_2d[i1[7],1:],yp2d[i1[7],1:]-yp2d[i1[7],0],'b-o',ms=4,label=r'$f_{k}$ in D-PANS')
axes[0].plot(fk_ref_2d[i1[0],1:],yp2d[i1[0],1:]-yp2d[i1[0],0],'r-o',ms=4,label=r'$f_{k} = C_{\mu}^{-0.5}(\frac{\Delta}{L_{t}})^{2/3}$')
axes[1].plot(fk_ref_2d[i1[1],1:],yp2d[i1[1],1:]-yp2d[i1[1],0],'r-o',ms=4,label=r'$f_{k} = C_{\mu}^{-0.5}(\frac{\Delta}{L_{t}})^{2/3}$')
axes[2].plot(fk_ref_2d[i1[2],1:],yp2d[i1[2],1:]-yp2d[i1[2],0],'r-o',ms=4,label=r'$f_{k} = C_{\mu}^{-0.5}(\frac{\Delta}{L_{t}})^{2/3}$')
axes[3].plot(fk_ref_2d[i1[3],1:],yp2d[i1[3],1:]-yp2d[i1[3],0],'r-o',ms=4,label=r'$f_{k} = C_{\mu}^{-0.5}(\frac{\Delta}{L_{t}})^{2/3}$')
axes[4].plot(fk_ref_2d[i1[4],1:],yp2d[i1[4],1:]-yp2d[i1[4],0],'r-o',ms=4,label=r'$f_{k} = C_{\mu}^{-0.5}(\frac{\Delta}{L_{t}})^{2/3}$')
axes[5].plot(fk_ref_2d[i1[5],1:],yp2d[i1[5],1:]-yp2d[i1[5],0],'r-o',ms=4,label=r'$f_{k} = C_{\mu}^{-0.5}(\frac{\Delta}{L_{t}})^{2/3}$')
axes[6].plot(fk_ref_2d[i1[6],1:],yp2d[i1[6],1:]-yp2d[i1[6],0],'r-o',ms=4,label=r'$f_{k} = C_{\mu}^{-0.5}(\frac{\Delta}{L_{t}})^{2/3}$')
axes[7].plot(fk_ref_2d[i1[7],1:],yp2d[i1[7],1:]-yp2d[i1[7],0],'r-o',ms=4,label=r'$f_{k} = C_{\mu}^{-0.5}(\frac{\Delta}{L_{t}})^{2/3}$')
for ax,xloc in zip(axes,xx):
    ax.plot(0.4*np.ones(20),np.linspace(0,0.2,20)-yp2d[i1[0],0],'k--',label=r'$f_{k}$ = 0.4')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(title='At X = '+str(xloc)+' m',title_fontsize=14,fontsize=14,loc='best')
    ax.tick_params(labelsize='medium')
    ax.grid(True, linestyle='-.')
    ax.set_ylim([0,0.2])
    ax.set_xlim([0,1])
    if (ax==axes[0] or ax==axes[4]):
        ax.set_ylabel("$y-y_{wall}$")
    if (ax==axes[4] or ax==axes[5] or ax==axes[6] or ax==axes[7]):
        ax.set_xlabel("$f_{k}$")
fig4.suptitle(r"$f_{k}$ computed at various locations",fontsize=22)

# %% Assignment T5 - Location of Interface in DES and DDES
# SA-DES model
C_DES_SA = 0.65
delta_SA = C_DES_SA * np.maximum.reduce([dx,dy,dz*np.ones([ni,nj])])
d_tilda = np.minimum(yp2d-yp2d[:,0][:,None],delta_SA)
switch_i_SA = (yp2d-yp2d[:,0][:,None]<=delta_SA).argmin(axis=1)
switch_wd_SA = np.zeros(len(switch_i_SA))
for i in range (ni):
    switch_wd_SA[i] = yp2d[i,switch_i_SA[i]] - yp2d[i,0]

# SST-DES model
C_DES_SST = 0.61
beta_star = 0.09    
omega_2d = (1/beta_star) * (diss_2d/k_model_2d)
Lt_SST = (1/beta_star) * ((k_model_2d**0.5)/omega_2d)
delta_SST = C_DES_SST * np.maximum.reduce([dx,dy,dz*np.ones([ni,nj])])
term1 = Lt_SST/delta_SST
F_SST = np.maximum(term1,1)
switch_i_SST = (term1<=1).argmin(axis=1)
switch_i_SST = np.where(switch_i_SST==0,109,switch_i_SST)
index_SST = np.where(switch_i_SST!=109)
switch_wd_SST = np.zeros(len(switch_i_SST))
switch_wd_SST2 = np.zeros(len(switch_i_SST))
for i in range (ni):
    switch_wd_SST[i] = yp2d[i,switch_i_SST[i]] - yp2d[i,0]
    switch_wd_SST2[i] = yp2d[i,switch_i_SST[i]]
switch_wd_SST1 = np.ma.masked_where(switch_wd_SST > 0.5, switch_wd_SST)
switch_wd_SST2 = np.ma.masked_where(switch_wd_SST2 > 0.5, switch_wd_SST2)

# DDES model
dist = yp2d[:]- yp2d[:,0][:,None]
term2 = 2*(Lt_SST[:,1:]/dist[:,1:])
term3 = (500*viscos)/(omega_2d[:,1:]*(dist[:,1:]**2))
eta = np.maximum(term2,term3)
F2 = np.tanh(eta**2)
term4 = term1[:,1:] * (1-F2)
switch_i_DDES = (term4<=1).argmin(axis=1)
switch_i_DDES = np.where(switch_i_DDES==0,109,switch_i_DDES)
index_DDES = np.where(switch_i_DDES!=109)
switch_wd_DDES = np.zeros(len(switch_i_SA))
switch_wd_DDES2 = np.zeros(len(switch_i_SA))
for i in range (ni):
    switch_wd_DDES[i] = yp2d[i,switch_i_DDES[i]] - yp2d[i,0]
    switch_wd_DDES2[i] = yp2d[i,switch_i_DDES[i]]
switch_wd_DDES1 = np.ma.masked_where(switch_wd_DDES > 0.5, switch_wd_DDES)
switch_wd_DDES2 = np.ma.masked_where(switch_wd_DDES2 > 0.5, switch_wd_DDES2)

# Length scales
ls_DES = C_DES_SA * np.maximum.reduce([dx,dy,dz*np.ones([ni,nj])])
ls_DPANS = (k_model_2d**1.5)/diss_2d

# finding boundary layer thickness
blt_i = np.zeros(ni)
blt_yc = np.zeros(ni)
blt_wd = np.zeros(ni)
for i in range(1,ni):
    blt_i[i] = np.max(np.where(((vis_2d[i,:]/viscos-1)<=1)==0))
    blt_yc[i] = yp2d[i,int(blt_i[i])]
    blt_wd[i] = blt_yc[i] - yp2d[i,0]

fig5,axes = plt.subplots(1,1,constrained_layout=True)
axes.plot(xp2d[1:-1,0],switch_wd,'b-o',ms=2,label='D-PANS model')
axes.plot(xp2d[1:-1,0],switch_wd_SA[1:-1],'r-o',ms=2,label='SA-DES model')
axes.plot(xp2d[1:-1,0],switch_wd_SST1[1:-1],'m-o',ms=2,label='SST-DES model')
axes.plot(xp2d[1:-1,0],switch_wd_DDES1[1:-1],'c-o',ms=2,label='DDES model')
axes.legend(fontsize=16,loc='best')
axes.xaxis.set_minor_locator(AutoMinorLocator())
axes.yaxis.set_minor_locator(AutoMinorLocator())
axes.tick_params(labelsize='medium')
axes.grid(True, linestyle='-.')
axes.set_xlabel('x (m)')
axes.set_ylabel(r'$y-y_{wall}$ (m)')
fig5.suptitle(r"Location of the switch (from RANS to LES)",fontsize=22)

fig5a,axes = plt.subplots(1,1,constrained_layout=True)
axes.plot(xp2d[1:-1,0],switch_wd+yp2d[1:-1,0],'b--',label='D-PANS model')
axes.plot(xp2d[1:-1,0],switch_wd_SA[1:-1]+yp2d[1:-1,0],'r--',label='SA-DES model')
axes.plot(xp2d[1:-1,0],switch_wd_SST2[1:-1],'m--',label='SST-DES model')
axes.plot(xp2d[1:-1,0],switch_wd_DDES2[1:-1],'c--',label='DDES model')
axes.plot(xp2d[:,0],yp2d[:,0],'k',lw=4,label='Lower boundary')
axes.legend(fontsize=16,loc='best')
axes.xaxis.set_minor_locator(AutoMinorLocator())
axes.yaxis.set_minor_locator(AutoMinorLocator())
axes.tick_params(labelsize='medium')
axes.grid(True, linestyle='-.')
axes.set_xlabel('x (m)')
axes.set_ylabel(r'$y$ (m)')
fig5a.suptitle(r"Location of the switch (from RANS to LES)",fontsize=22)

fig6,axes = plt.subplots(2,4,constrained_layout=True)
axes = axes.reshape(-1)
axes[0].plot(ls_DES[i1[1],1:],yp2d[i1[1],1:]-yp2d[i1[1],0],'b-o',ms=4,label=r'$l_{DES}$')
axes[1].plot(ls_DES[i1[2],1:],yp2d[i1[2],1:]-yp2d[i1[2],0],'b-o',ms=4,label=r'$l_{DES}$')
axes[2].plot(ls_DES[i1[3],1:],yp2d[i1[3],1:]-yp2d[i1[3],0],'b-o',ms=4,label=r'$l_{DES}$')
axes[3].plot(ls_DES[i1[4],1:],yp2d[i1[4],1:]-yp2d[i1[4],0],'b-o',ms=4,label=r'$l_{DES}$')
axes[4].plot(ls_DES[i1[5],1:],yp2d[i1[5],1:]-yp2d[i1[5],0],'b-o',ms=4,label=r'$l_{DES}$')
axes[5].plot(ls_DES[i1[6],1:],yp2d[i1[6],1:]-yp2d[i1[6],0],'b-o',ms=4,label=r'$l_{DES}$')
axes[6].plot(ls_DES[i1[7],1:],yp2d[i1[7],1:]-yp2d[i1[7],0],'b-o',ms=4,label=r'$l_{DES}$')
axes[7].plot(ls_DES[i1[8],1:],yp2d[i1[8],1:]-yp2d[i1[8],0],'b-o',ms=4,label=r'$l_{DES}$')
axes[0].plot(ls_DPANS[i1[1],1:],yp2d[i1[1],1:]-yp2d[i1[1],0],'r-o',ms=4,label=r'$l_{D-PANS}$')
axes[1].plot(ls_DPANS[i1[2],1:],yp2d[i1[2],1:]-yp2d[i1[2],0],'r-o',ms=4,label=r'$l_{D-PANS}$')
axes[2].plot(ls_DPANS[i1[3],1:],yp2d[i1[3],1:]-yp2d[i1[3],0],'r-o',ms=4,label=r'$l_{D-PANS}$')
axes[3].plot(ls_DPANS[i1[4],1:],yp2d[i1[4],1:]-yp2d[i1[4],0],'r-o',ms=4,label=r'$l_{D-PANS}$')
axes[4].plot(ls_DPANS[i1[5],1:],yp2d[i1[5],1:]-yp2d[i1[5],0],'r-o',ms=4,label=r'$l_{D-PANS}$')
axes[5].plot(ls_DPANS[i1[6],1:],yp2d[i1[6],1:]-yp2d[i1[6],0],'r-o',ms=4,label=r'$l_{D-PANS}$')
axes[6].plot(ls_DPANS[i1[7],1:],yp2d[i1[7],1:]-yp2d[i1[7],0],'r-o',ms=4,label=r'$l_{D-PANS}$')
axes[7].plot(ls_DPANS[i1[8],1:],yp2d[i1[8],1:]-yp2d[i1[8],0],'r-o',ms=4,label=r'$l_{D-PANS}$')
axes[0].plot(np.linspace(np.min(ls_DPANS[i1[1],1:]),np.max(ls_DES[i1[1],1:]),10),blt_wd[i1[1]]*np.ones(10),'k--',ms=4,label='BL thickness')
axes[1].plot(np.linspace(np.min(ls_DPANS[i1[2],1:]),np.max(ls_DES[i1[2],1:]),10),blt_wd[i1[2]]*np.ones(10),'k--',ms=4,label='BL thickness')
axes[2].plot(np.linspace(np.min(ls_DPANS[i1[3],1:]),np.max(ls_DES[i1[3],1:]),10),blt_wd[i1[3]]*np.ones(10),'k--',ms=4,label='BL thickness')
axes[3].plot(np.linspace(np.min(ls_DPANS[i1[4],1:]),np.max(ls_DES[i1[4],1:]),10),blt_wd[i1[4]]*np.ones(10),'k--',ms=4,label='BL thickness')
axes[4].plot(np.linspace(np.min(ls_DPANS[i1[5],1:]),np.max(ls_DES[i1[5],1:]),10),blt_wd[i1[5]]*np.ones(10),'k--',ms=4,label='BL thickness')
axes[5].plot(np.linspace(np.min(ls_DPANS[i1[6],1:]),np.max(ls_DES[i1[6],1:]),10),blt_wd[i1[6]]*np.ones(10),'k--',ms=4,label='BL thickness')
axes[6].plot(np.linspace(np.min(ls_DPANS[i1[7],1:]),np.max(ls_DES[i1[7],1:]),10),blt_wd[i1[7]]*np.ones(10),'k--',ms=4,label='BL thickness')
axes[7].plot(np.linspace(np.min(ls_DPANS[i1[8],1:]),np.max(ls_DES[i1[8],1:]),10),blt_wd[i1[8]]*np.ones(10),'k--',ms=4,label='BL thickness')
for ax,xloc in zip(axes,xx[1:]):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(title='At X = '+str(xloc)+' m',title_fontsize=14,fontsize=14,loc='best')
    ax.tick_params(labelsize='medium')
    ax.grid(True, linestyle='-.')
    ax.set_xlim([0,0.015])
    ax.set_ylim([0,0.2])
    if ax==axes[7]:
        ax.set_xlim([0,0.03])
        ax.set_ylim([0,0.3])
    if (ax==axes[0] or ax==axes[4]):
        ax.set_ylabel("$y-y_{wall}$")
    if (ax==axes[4] or ax==axes[5] or ax==axes[6] or ax==axes[7]):
        ax.set_xlabel("Length scale")
fig6.suptitle(r"Length scales computed at various locations",fontsize=22)