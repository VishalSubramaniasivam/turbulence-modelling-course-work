import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from dphidx_dy import dphidx_dy
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
plt.rcParams.update({'font.size': 22})

re =9.36e+5
viscos =1/re

xy_hump_fine = np.loadtxt("xy_hump.dat")
x=xy_hump_fine[:,0]
y=xy_hump_fine[:,1]

ni=314
nj=122

nim1=ni-1
njm1=nj-1

# Read Data File (from paper)
vectz=np.genfromtxt("vectz_aiaa_paper.dat",comments="%")
ntstep=vectz[0]
n=len(vectz)

# Read Data File (from journal)
vectz_aiaa = np.genfromtxt("vectz_aiaa_journal.dat",comments="%")

#            write(48,*)uvec(i,j)
#            write(48,*)vvec(i,j)
#            write(48,*)dummy(i,j)
#            write(48,*)uvec2(i,j)
#            write(48,*)vvec2(i,j)
#            write(48,*)wvec2(i,j)
#            write(48,*)uvvec(i,j)
#            write(48,*)p2D(i,j)
#            write(48,*)rk2D(i,j)
#            write(48,*)vis2D(i,j)  
#            write(48,*)dissp2D(i,j)
#            write(48,*)uvturb(i,j)



nn=12
nst=0
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
iuv_model=range(nst+12,n,nn)

# From paper
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
uv_model=vectz[iuv_model]/ntstep
# From journal
u1=vectz_aiaa[iu]/ntstep
v1=vectz_aiaa[iv]/ntstep
fk1=vectz_aiaa[ifk]/ntstep
uu1=vectz_aiaa[iuu]/ntstep
vv1=vectz_aiaa[ivv]/ntstep
ww1=vectz_aiaa[iww]/ntstep
uv1=vectz_aiaa[iuv]/ntstep
p1=vectz_aiaa[ip]/ntstep
k_model1=vectz_aiaa[ik]/ntstep
vis1=vectz_aiaa[ivis]/ntstep
diss1=vectz_aiaa[idiss]/ntstep
uv1_model=vectz_aiaa[iuv_model]/ntstep


# uu is total inst. velocity squared. Hence the resolved turbulent resolved stresses are obtained as
# From paper
uu=uu-u**2
vv=vv-v**2
uv=uv-u*v
# From journal
uu1=uu1-u**2
vv1=vv1-v**2
uv1=uv1-u*v

# From paper
p_2d=np.reshape(p,(ni,nj))
u_2d=np.reshape(u,(ni,nj))
v_2d=np.reshape(v,(ni,nj))
fk_2d=np.reshape(fk,(ni,nj))
uu_2d=np.reshape(uu,(ni,nj))
uv_2d=np.reshape(uv,(ni,nj))
vv_2d=np.reshape(vv,(ni,nj))
ww_2d=np.reshape(ww,(ni,nj))
uv_model_2d=np.reshape(uv_model,(ni,nj))
k_model_2d=np.reshape(k_model,(ni,nj))
vis_2d=np.reshape(vis,(ni,nj)) #this is to total viscosity, i.e. vis_tot=vis+vis_turb
diss_2d=np.reshape(diss,(ni,nj)) #this is to total viscosity, i.e. vis_tot=vis+vis_turb
x_2d=np.transpose(np.reshape(x,(nj,ni)))
y_2d=np.transpose(np.reshape(y,(nj,ni)))
# From journal
p1_2d=np.reshape(p1,(ni,nj))
u1_2d=np.reshape(u1,(ni,nj))
v1_2d=np.reshape(v1,(ni,nj))
fk1_2d=np.reshape(fk1,(ni,nj))
uu1_2d=np.reshape(uu1,(ni,nj))
uv1_2d=np.reshape(uv1,(ni,nj))
vv1_2d=np.reshape(vv1,(ni,nj))
ww1_2d=np.reshape(ww1,(ni,nj))
uv1_model_2d=np.reshape(uv1_model,(ni,nj))
k_model1_2d=np.reshape(k_model1,(ni,nj))
vis1_2d=np.reshape(vis1,(ni,nj)) #this is to total viscosity, i.e. vis_tot=vis+vis_turb
diss1_2d=np.reshape(diss1,(ni,nj)) #this is to total viscosity, i.e. vis_tot=vis+vis_turb

# set fk_2d=1 at upper boundary
fk_2d[:,nj-1]=fk_2d[:,nj-2]
fk1_2d[:,nj-1]=fk_2d[:,nj-2]

x065_off=np.genfromtxt("x065_off.dat",comments="%")
x066_off=np.genfromtxt("x066_off.dat",comments="%")
x080_off=np.genfromtxt("x080_off.dat",comments="%")
x090_off=np.genfromtxt("x090_off.dat",comments="%")
x100_off=np.genfromtxt("x100_off.dat",comments="%")
x110_off=np.genfromtxt("x110_off.dat",comments="%")
x120_off=np.genfromtxt("x120_off.dat",comments="%")
x130_off=np.genfromtxt("x130_off.dat",comments="%")

# face cordinates for dphidx_dy function
#xf2d = np.zeros((ni,nj))
#yf2d = np.zeros((ni,nj))
xf2d = 0.5 * (x_2d[1:-2,:] + x_2d[2:-1,:])
yf2d = 0.5 * (y_2d[:,1:-2] + y_2d[:,2:-1])
#xf2d[] = x_2d[]
#yf2d[] = y_2d[]

# the funtion dphidx_dy wants x and y arrays to be one cell smaller than u2d. Hence I take away the last row and column below
x_2d_new=np.delete(x_2d,-1,0)
x_2d_new=np.delete(x_2d_new,-1,1)
y_2d_new=np.delete(y_2d,-1,0)
y_2d_new=np.delete(y_2d_new,-1,1)
# compute the gradient
dudx,dudy=dphidx_dy(x_2d_new,y_2d_new,u_2d)
dvdx,dvdy=dphidx_dy(x_2d_new,y_2d_new,v_2d)

#*************************
# Mesh for visualization
a = np.transpose(x_2d)
b = np.transpose(y_2d)
fig16 = plt.figure()
axes = fig16.add_axes([0.1,0.1,0.8,0.8])
axes.plot(x_2d,y_2d,color='black')
axes.plot(a,b,color='black')
axes.tick_params(axis="x", labelsize=15)
axes.tick_params(axis="y", labelsize=15)
axes.set_xlabel('x [m]',fontsize=20)
axes.set_ylabel('y [m]',fontsize=20)
axes.set_title('Computational mesh in the domain',fontsize=24)

#*************************
# plot u
fig1 = plt.figure("Figure 1")
plt.clf() #clear the figure
xx=0.65;
i1 = (np.abs(xx-x_2d[:,1])).argmin()  # find index which closest fits xx
plt.plot(u_2d[i1,:],y_2d[i1,:]-y_2d[i1,0],'b-')
plt.plot(x065_off[:,2],x065_off[:,1]-y_2d[i1,0],'bo')
plt.xlabel("$U$")
plt.ylabel("$y-y_{wall}$")
plt.title("$x=0.65$")
plt.axis([0, 1.3,0,0.2])
#plt.savefig('u065_hump_python.eps')

#*************************
# plot vv
fig2 = plt.figure("Figure 2")
plt.clf() #clear the figure
xx=0.65;
i1 = (np.abs(xx-x_2d[:,1])).argmin()  # find index which closest fits xx
plt.plot(vv_2d[i1,:],y_2d[i1,:]-y_2d[i1,0],'b-')
plt.plot(x065_off[:,5],x065_off[:,1]-y_2d[i1,0],'bo')
plt.xlabel("$\overline{v'v'}$")
plt.ylabel("$y-y_{wall}$")
plt.title("$x=0.65$")
plt.axis([0, 0.01,0,0.2])
#plt.savefig('vv065_hump_python.eps')

#*************************
# plot uu
fig3 = plt.figure("Figure 3")
plt.clf() #clear the figure
xx=0.65;
i1 = (np.abs(xx-x_2d[:,1])).argmin()  # find index which closest fits xx
plt.plot(uu_2d[i1,:],y_2d[i1,:]-y_2d[i1,0],'b-')
plt.plot(x065_off[:,4],x065_off[:,1]-y_2d[i1,0],'bo')
plt.xlabel("$\overline{u'u'}$")
plt.ylabel("$y-y_{wall}$")
plt.title("$x=0.65$")
plt.axis([0, 0.01,0,0.2])
#plt.savefig('uu065_hump_python.eps')

################################ contoue plot
fig4 = plt.figure("Figure 4")
plt.clf() #clear the figure
plt.contourf(x_2d,y_2d,uu_2d, 50)
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.clim(0,0.05)
plt.axis([0.6,1.5,0,1])
plt.title("contour $\overline{u'u'}$")
#plt.savefig('piso_python.eps')

################################ vector plot
fig5 = plt.figure("Figure 5")
plt.clf() #clear the figure
k=6# plot every forth vector
ss=3.2 #vector length
plt.quiver(x_2d[::k,::k],y_2d[::k,::k],u_2d[::k,::k],v_2d[::k,::k],width=0.01)
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.axis([0.6,1.5,0,1])
plt.title("vector plot")
#plt.savefig('vect_python.eps')

# %% Assignment T1 - Discretization Schemes
# Paper Data [203]
## Plot vv
xx = np.array([0,0.65,0.8,0.9,1,1.1,1.2,1.3])
i1 = list(map(lambda m : (np.abs(m-x_2d[:,1])).argmin(),xx))
fig6,axes = plt.subplots(2,4,constrained_layout=True)
axes = axes.reshape(-1)
axes[0].plot(vv_2d[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'b-o',ms=4,label='x = 0 m')
axes[1].plot(vv_2d[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'g-v',ms=4,label='x = 0.65 m')
axes[2].plot(vv_2d[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'r-1',ms=4,label='x = 0.8 m')
axes[3].plot(vv_2d[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'c-s',ms=4,label='x = 0.9 m')
axes[4].plot(vv_2d[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'m-+',ms=4,label='x = 1 m')
axes[5].plot(vv_2d[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'y-*',ms=4,label='x = 1.1 m')
axes[6].plot(vv_2d[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'k-d',ms=4,label='x = 1.2 m')
axes[7].plot(vv_2d[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],label='x = 1.3 m')
for ax in axes:
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(fontsize=14,loc='best')
    ax.tick_params(labelsize='medium')
    ax.grid(True, linestyle='-.')
    ax.set_ylabel("$y-y_{wall}$")
    if (ax==axes[4] or ax==axes[5] or ax==axes[6] or ax==axes[7]):
        ax.set_xlabel("$\overline{v'v'}$")
fig6.suptitle("$<\overline{v'v'}>$ at various locations",fontsize=22)

## Plot uv
fig7,axes = plt.subplots(2,4,constrained_layout=True)
axes = axes.reshape(-1)
axes[0].plot(uv_2d[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'b-o',ms=4,label='x = 0 m')
axes[1].plot(uv_2d[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'g-v',ms=4,label='x = 0.65 m')
axes[2].plot(uv_2d[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'r-1',ms=4,label='x = 0.8 m')
axes[3].plot(uv_2d[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'c-s',ms=4,label='x = 0.9 m')
axes[4].plot(uv_2d[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'m-+',ms=4,label='x = 1 m')
axes[5].plot(uv_2d[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'y-*',ms=4,label='x = 1.1 m')
axes[6].plot(uv_2d[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'k-d',ms=4,label='x = 1.2 m')
axes[7].plot(uv_2d[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],label='x = 1.3 m')
for ax in axes:
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(fontsize=14,loc='best')
    ax.tick_params(labelsize='medium')
    ax.grid(True, linestyle='-.')
    ax.set_ylabel("$y-y_{wall}$")
    if (ax==axes[4] or ax==axes[5] or ax==axes[6] or ax==axes[7]):
        ax.set_xlabel("$\overline{u'v'}$")
fig7.suptitle("<$\overline{u'v'}$> at various locations",fontsize=22)

# Journal Data [169]
## Plot vv
fig8,axes = plt.subplots(2,4,constrained_layout=True)
axes = axes.reshape(-1)
axes[0].plot(vv1_2d[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'b-o',ms=4,label='x = 0 m')
axes[1].plot(vv1_2d[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'g-v',ms=4,label='x = 0.65 m')
axes[2].plot(vv1_2d[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'r-1',ms=4,label='x = 0.8 m')
axes[3].plot(vv1_2d[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'c-s',ms=4,label='x = 0.9 m')
axes[4].plot(vv1_2d[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'m-+',ms=4,label='x = 1 m')
axes[5].plot(vv1_2d[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'y-*',ms=4,label='x = 1.1 m')
axes[6].plot(vv1_2d[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'k-d',ms=4,label='x = 1.2 m')
axes[7].plot(vv1_2d[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],label='x = 1.3 m')
for ax in axes:
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(fontsize=14,loc='best')
    ax.tick_params(labelsize='medium')
    ax.grid(True, linestyle='-.')
    ax.set_ylabel("$y-y_{wall}$")
    if (ax==axes[4] or ax==axes[5] or ax==axes[6] or ax==axes[7]):
        ax.set_xlabel("$\overline{v'v'}$")
fig8.suptitle("$<\overline{v'v'}>$ at various locations",fontsize=22)

## Plot uv
fig9,axes = plt.subplots(2,4,constrained_layout=True)
axes = axes.reshape(-1)
axes[0].plot(uv1_2d[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'b-o',ms=4,label='x = 0 m')
axes[1].plot(uv1_2d[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'g-v',ms=4,label='x = 0.65 m')
axes[2].plot(uv1_2d[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'r-1',ms=4,label='x = 0.8 m')
axes[3].plot(uv1_2d[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'c-s',ms=4,label='x = 0.9 m')
axes[4].plot(uv1_2d[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'m-+',ms=4,label='x = 1 m')
axes[5].plot(uv1_2d[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'y-*',ms=4,label='x = 1.1 m')
axes[6].plot(uv1_2d[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'k-d',ms=4,label='x = 1.2 m')
axes[7].plot(uv1_2d[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],label='x = 1.3 m')
for ax in axes:
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(fontsize=14,loc='best')
    ax.tick_params(labelsize='medium')
    ax.grid(True, linestyle='-.')
    ax.set_ylabel("$y-y_{wall}$")
    if (ax==axes[4] or ax==axes[5] or ax==axes[6] or ax==axes[7]):
        ax.set_xlabel("$\overline{u'v'}$")
fig9.suptitle("<$\overline{v'v'}$> at various locations",fontsize=22)

# Comparision of both dataset
## Plot uu
fig10,axes = plt.subplots(2,4,constrained_layout=True)
axes = axes.reshape(-1)
axes[0].plot(uu_2d[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'r-o',ms=4,label='Paper Data [203]')
axes[1].plot(uu_2d[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'r-o',ms=4,label='Paper Data [203]')
axes[2].plot(uu_2d[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'r-o',ms=4,label='Paper Data [203]')
axes[3].plot(uu_2d[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'r-o',ms=4,label='Paper Data [203]')
axes[4].plot(uu_2d[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'r-o',ms=4,label='Paper Data [203]')
axes[5].plot(uu_2d[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'r-o',ms=4,label='Paper Data [203]')
axes[6].plot(uu_2d[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'r-o',ms=4,label='Paper Data [203]')
axes[7].plot(uu_2d[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'r-o',ms=4,label='Paper Data [203]')
axes[0].plot(uu1_2d[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'b-s',ms=4,label='Journal Data [169]')
axes[1].plot(uu1_2d[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'b-s',ms=4,label='Journal Data [169]')
axes[2].plot(uu1_2d[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'b-s',ms=4,label='Journal Data [169]')
axes[3].plot(uu1_2d[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'b-s',ms=4,label='Journal Data [169]')
axes[4].plot(uu1_2d[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'b-s',ms=4,label='Journal Data [169]')
axes[5].plot(uu1_2d[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'b-s',ms=4,label='Journal Data [169]')
axes[6].plot(uu1_2d[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'b-s',ms=4,label='Journal Data [169]')
axes[7].plot(uu1_2d[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'b-s',ms=4,label='Journal Data [169]')
axes[1].plot(x065_off[:,4],x065_off[:,1]-y_2d[i1[1],0],'k-*',ms=4,label='Experimental Data')
axes[2].plot(x080_off[:,4],x080_off[:,1]-y_2d[i1[2],0],'k-*',ms=4,label='Experimental Data')
axes[3].plot(x090_off[:,4],x090_off[:,1]-y_2d[i1[3],0],'k-*',ms=4,label='Experimental Data')
axes[4].plot(x100_off[:,4],x100_off[:,1]-y_2d[i1[4],0],'k-*',ms=4,label='Experimental Data')
axes[5].plot(x110_off[:,4],x110_off[:,1]-y_2d[i1[5],0],'k-*',ms=4,label='Experimental Data')
axes[6].plot(x120_off[:,4],x120_off[:,1]-y_2d[i1[6],0],'k-*',ms=4,label='Experimental Data')
axes[7].plot(x130_off[:,4],x130_off[:,1]-y_2d[i1[7],0],'k-*',ms=4,label='Experimental Data')
for ax,xloc in zip(axes,xx):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(title='At X = '+str(xloc)+' m',title_fontsize=14,fontsize=14,loc='best')
    ax.tick_params(labelsize='medium')
    ax.grid(True, linestyle='-.')
    ax.set_ylim([0,0.25])
    if (ax==axes[4] or ax==axes[5] or ax==axes[6] or ax==axes[7]):
        ax.set_xlabel("$<\overline{u'u'}>$")
    if (ax==axes[0] or ax==axes[4]):
        ax.set_ylabel("$y-y_{wall}$")
fig10.suptitle("$<\overline{u'u'}>$ at various locations",fontsize=22)

## Plot vv
fig11,axes = plt.subplots(2,4,constrained_layout=True)
axes = axes.reshape(-1)
axes[0].plot(vv_2d[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'r-o',ms=4,label='Paper Data [203]')
axes[1].plot(vv_2d[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'r-o',ms=4,label='Paper Data [203]')
axes[2].plot(vv_2d[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'r-o',ms=4,label='Paper Data [203]')
axes[3].plot(vv_2d[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'r-o',ms=4,label='Paper Data [203]')
axes[4].plot(vv_2d[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'r-o',ms=4,label='Paper Data [203]')
axes[5].plot(vv_2d[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'r-o',ms=4,label='Paper Data [203]')
axes[6].plot(vv_2d[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'r-o',ms=4,label='Paper Data [203]')
axes[7].plot(vv_2d[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'r-o',ms=4,label='Paper Data [203]')
axes[0].plot(vv1_2d[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'b-s',ms=4,label='Journal Data [169]')
axes[1].plot(vv1_2d[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'b-s',ms=4,label='Journal Data [169]')
axes[2].plot(vv1_2d[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'b-s',ms=4,label='Journal Data [169]')
axes[3].plot(vv1_2d[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'b-s',ms=4,label='Journal Data [169]')
axes[4].plot(vv1_2d[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'b-s',ms=4,label='Journal Data [169]')
axes[5].plot(vv1_2d[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'b-s',ms=4,label='Journal Data [169]')
axes[6].plot(vv1_2d[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'b-s',ms=4,label='Journal Data [169]')
axes[7].plot(vv1_2d[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'b-s',ms=4,label='Journal Data [169]')
axes[1].plot(x065_off[:,5],x065_off[:,1]-y_2d[i1[1],0],'k-*',ms=4,label='Experimental Data')
axes[2].plot(x080_off[:,5],x080_off[:,1]-y_2d[i1[2],0],'k-*',ms=4,label='Experimental Data')
axes[3].plot(x090_off[:,5],x090_off[:,1]-y_2d[i1[3],0],'k-*',ms=4,label='Experimental Data')
axes[4].plot(x100_off[:,5],x100_off[:,1]-y_2d[i1[4],0],'k-*',ms=4,label='Experimental Data')
axes[5].plot(x110_off[:,5],x110_off[:,1]-y_2d[i1[5],0],'k-*',ms=4,label='Experimental Data')
axes[6].plot(x120_off[:,5],x120_off[:,1]-y_2d[i1[6],0],'k-*',ms=4,label='Experimental Data')
axes[7].plot(x130_off[:,5],x130_off[:,1]-y_2d[i1[7],0],'k-*',ms=4,label='Experimental Data')
for ax,xloc in zip(axes,xx):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(title='At X = '+str(xloc)+' m',title_fontsize=14,fontsize=14,loc='best')
    ax.tick_params(labelsize='medium')
    ax.grid(True, linestyle='-.')
    ax.set_ylim([0,0.25])
    if (ax==axes[4] or ax==axes[5] or ax==axes[6] or ax==axes[7]):
        ax.set_xlabel("$<\overline{v'v'}>$")
    if (ax==axes[0] or ax==axes[4]):
        ax.set_ylabel("$y-y_{wall}$")
fig11.suptitle("$<\overline{v'v'}>$ at various locations",fontsize=22)

## Plot uv
fig12,axes = plt.subplots(2,4,constrained_layout=True)
axes = axes.reshape(-1)
axes[0].plot(uv_2d[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'r-o',ms=4,label='Paper Data [203]')
axes[1].plot(uv_2d[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'r-o',ms=4,label='Paper Data [203]')
axes[2].plot(uv_2d[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'r-o',ms=4,label='Paper Data [203]')
axes[3].plot(uv_2d[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'r-o',ms=4,label='Paper Data [203]')
axes[4].plot(uv_2d[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'r-o',ms=4,label='Paper Data [203]')
axes[5].plot(uv_2d[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'r-o',ms=4,label='Paper Data [203]')
axes[6].plot(uv_2d[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'r-o',ms=4,label='Paper Data [203]')
axes[7].plot(uv_2d[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'r-o',ms=4,label='Paper Data [203]')
axes[0].plot(uv1_2d[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'b-s',ms=4,label='Journal Data [169]')
axes[1].plot(uv1_2d[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'b-s',ms=4,label='Journal Data [169]')
axes[2].plot(uv1_2d[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'b-s',ms=4,label='Journal Data [169]')
axes[3].plot(uv1_2d[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'b-s',ms=4,label='Journal Data [169]')
axes[4].plot(uv1_2d[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'b-s',ms=4,label='Journal Data [169]')
axes[5].plot(uv1_2d[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'b-s',ms=4,label='Journal Data [169]')
axes[6].plot(uv1_2d[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'b-s',ms=4,label='Journal Data [169]')
axes[7].plot(uv1_2d[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'b-s',ms=4,label='Journal Data [169]')
axes[1].plot(x065_off[:,6],x065_off[:,1]-y_2d[i1[1],0],'k-*',ms=4,label='Experimental Data')
axes[2].plot(x080_off[:,6],x080_off[:,1]-y_2d[i1[2],0],'k-*',ms=4,label='Experimental Data')
axes[3].plot(x090_off[:,6],x090_off[:,1]-y_2d[i1[3],0],'k-*',ms=4,label='Experimental Data')
axes[4].plot(x100_off[:,6],x100_off[:,1]-y_2d[i1[4],0],'k-*',ms=4,label='Experimental Data')
axes[5].plot(x110_off[:,6],x110_off[:,1]-y_2d[i1[5],0],'k-*',ms=4,label='Experimental Data')
axes[6].plot(x120_off[:,6],x120_off[:,1]-y_2d[i1[6],0],'k-*',ms=4,label='Experimental Data')
axes[7].plot(x130_off[:,6],x130_off[:,1]-y_2d[i1[7],0],'k-*',ms=4,label='Experimental Data')
for ax,xloc in zip(axes,xx):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(title='At X = '+str(xloc)+' m',title_fontsize=14,fontsize=14,loc='best')
    ax.tick_params(labelsize='medium')
    ax.grid(True, linestyle='-.')
    ax.set_ylim([0,0.2])
    if (ax==axes[4] or ax==axes[5] or ax==axes[6] or ax==axes[7]):
        ax.set_xlabel("$<\overline{u'v'}>$")
    if (ax==axes[0] or ax==axes[4]):
        ax.set_ylabel("$y-y_{wall}$")
fig12.suptitle("$<\overline{u'v'}>$ at various locations",fontsize=22)

## Plot u
fig13,axes = plt.subplots(2,4,constrained_layout=True)
axes = axes.reshape(-1)
axes[0].plot(u_2d[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'r-o',ms=4,label='Paper Data [203]')
axes[1].plot(u_2d[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'r-o',ms=4,label='Paper Data [203]')
axes[2].plot(u_2d[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'r-o',ms=4,label='Paper Data [203]')
axes[3].plot(u_2d[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'r-o',ms=4,label='Paper Data [203]')
axes[4].plot(u_2d[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'r-o',ms=4,label='Paper Data [203]')
axes[5].plot(u_2d[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'r-o',ms=4,label='Paper Data [203]')
axes[6].plot(u_2d[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'r-o',ms=4,label='Paper Data [203]')
axes[7].plot(u_2d[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'r-o',ms=4,label='Paper Data [203]')
axes[0].plot(u1_2d[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'b-s',ms=4,label='Journal Data [169]')
axes[1].plot(u1_2d[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'b-s',ms=4,label='Journal Data [169]')
axes[2].plot(u1_2d[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'b-s',ms=4,label='Journal Data [169]')
axes[3].plot(u1_2d[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'b-s',ms=4,label='Journal Data [169]')
axes[4].plot(u1_2d[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'b-s',ms=4,label='Journal Data [169]')
axes[5].plot(u1_2d[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'b-s',ms=4,label='Journal Data [169]')
axes[6].plot(u1_2d[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'b-s',ms=4,label='Journal Data [169]')
axes[7].plot(u1_2d[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'b-s',ms=4,label='Journal Data [169]')
axes[1].plot(x065_off[:,2],x065_off[:,1]-y_2d[i1[1],0],'k-*',ms=4,label='Experimental Data')
axes[2].plot(x080_off[:,2],x080_off[:,1]-y_2d[i1[2],0],'k-*',ms=4,label='Experimental Data')
axes[3].plot(x090_off[:,2],x090_off[:,1]-y_2d[i1[3],0],'k-*',ms=4,label='Experimental Data')
axes[4].plot(x100_off[:,2],x100_off[:,1]-y_2d[i1[4],0],'k-*',ms=4,label='Experimental Data')
axes[5].plot(x110_off[:,2],x110_off[:,1]-y_2d[i1[5],0],'k-*',ms=4,label='Experimental Data')
axes[6].plot(x120_off[:,2],x120_off[:,1]-y_2d[i1[6],0],'k-*',ms=4,label='Experimental Data')
axes[7].plot(x130_off[:,2],x130_off[:,1]-y_2d[i1[7],0],'k-*',ms=4,label='Experimental Data')
for ax,xloc in zip(axes,xx):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(title='At X = '+str(xloc)+' m',title_fontsize=14,fontsize=14,loc='best')
    ax.tick_params(labelsize='medium')
    ax.grid(True, linestyle='-.')
    ax.set_ylim([0,0.15])
    if (ax==axes[4] or ax==axes[5] or ax==axes[6] or ax==axes[7]):
        ax.set_xlabel("$<\overline{u}>$")
    if (ax==axes[0] or ax==axes[4]):
        ax.set_ylabel("$y-y_{wall}$")
fig13.suptitle("$<\overline{u}>$ at various locations",fontsize=22)

## Plot v
fig14,axes = plt.subplots(2,4,constrained_layout=True)
axes = axes.reshape(-1)
axes[0].plot(v_2d[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'r-o',ms=4,label='Paper Data [203]')
axes[1].plot(v_2d[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'r-o',ms=4,label='Paper Data [203]')
axes[2].plot(v_2d[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'r-o',ms=4,label='Paper Data [203]')
axes[3].plot(v_2d[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'r-o',ms=4,label='Paper Data [203]')
axes[4].plot(v_2d[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'r-o',ms=4,label='Paper Data [203]')
axes[5].plot(v_2d[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'r-o',ms=4,label='Paper Data [203]')
axes[6].plot(v_2d[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'r-o',ms=4,label='Paper Data [203]')
axes[7].plot(v_2d[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'r-o',ms=4,label='Paper Data [203]')
axes[0].plot(v1_2d[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'b-s',ms=4,label='Journal Data [169]')
axes[1].plot(v1_2d[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'b-s',ms=4,label='Journal Data [169]')
axes[2].plot(v1_2d[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'b-s',ms=4,label='Journal Data [169]')
axes[3].plot(v1_2d[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'b-s',ms=4,label='Journal Data [169]')
axes[4].plot(v1_2d[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'b-s',ms=4,label='Journal Data [169]')
axes[5].plot(v1_2d[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'b-s',ms=4,label='Journal Data [169]')
axes[6].plot(v1_2d[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'b-s',ms=4,label='Journal Data [169]')
axes[7].plot(v1_2d[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'b-s',ms=4,label='Journal Data [169]')
axes[1].plot(x065_off[:,3],x065_off[:,1]-y_2d[i1[1],0],'k-*',ms=4,label='Experimental Data')
axes[2].plot(x080_off[:,3],x080_off[:,1]-y_2d[i1[2],0],'k-*',ms=4,label='Experimental Data')
axes[3].plot(x090_off[:,3],x090_off[:,1]-y_2d[i1[3],0],'k-*',ms=4,label='Experimental Data')
axes[4].plot(x100_off[:,3],x100_off[:,1]-y_2d[i1[4],0],'k-*',ms=4,label='Experimental Data')
axes[5].plot(x110_off[:,3],x110_off[:,1]-y_2d[i1[5],0],'k-*',ms=4,label='Experimental Data')
axes[6].plot(x120_off[:,3],x120_off[:,1]-y_2d[i1[6],0],'k-*',ms=4,label='Experimental Data')
axes[7].plot(x130_off[:,3],x130_off[:,1]-y_2d[i1[7],0],'k-*',ms=4,label='Experimental Data')
for ax,xloc in zip(axes,xx):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(title='At X = '+str(xloc)+' m',title_fontsize=14,fontsize=14,loc='best')
    ax.tick_params(labelsize='medium')
    ax.grid(True, linestyle='-.')
    ax.set_ylim([0,0.15])
    if (ax==axes[4] or ax==axes[5] or ax==axes[6] or ax==axes[7]):
        ax.set_xlabel("$<\overline{v}>$")
    if (ax==axes[0] or ax==axes[4]):
        ax.set_ylabel("$y-y_{wall}$")
fig14.suptitle("$<\overline{v}>$ at various locations",fontsize=22)

# %% Assignment T2 - Modeled Turbulent Shear Stress
visturb_2d = vis_2d - viscos
visturb1_2d = vis1_2d - viscos
uv_model_2d_calc = -visturb_2d * (dudy+dvdx)

# Finding boundary layer thickness at the measuring points (Paper Data)
tbl065 = np.max(np.where(((vis_2d[i1[1],:]/viscos-1)<=1) == 0))
tbl080 = np.max(np.where(((vis_2d[i1[2],:]/viscos-1)<=1) == 0))
tbl090 = np.max(np.where(((vis_2d[i1[3],:]/viscos-1)<=1) == 0))
tbl100 = np.max(np.where(((vis_2d[i1[4],:]/viscos-1)<=1) == 0))
tbl110 = np.max(np.where(((vis_2d[i1[5],:]/viscos-1)<=1) == 0))
tbl120 = np.max(np.where(((vis_2d[i1[6],:]/viscos-1)<=1) == 0))
tbl130 = np.max(np.where(((vis_2d[i1[7],:]/viscos-1)<=1) == 0))

blt060 = y_2d[i1[1],tbl065]
blt080 = y_2d[i1[2],tbl080]
blt090 = y_2d[i1[3],tbl090]
blt100 = y_2d[i1[4],tbl100]
blt110 = y_2d[i1[5],tbl110]
blt120 = y_2d[i1[6],tbl120]
blt130 = y_2d[i1[7],tbl130]

# Finding boundary layer thickness at the measuring points (Journal Data)
tbl1065 = np.max(np.where(((vis1_2d[i1[1],:]/viscos-1)<=1) == 0))
tbl1080 = np.max(np.where(((vis1_2d[i1[2],:]/viscos-1)<=1) == 0))
tbl1090 = np.max(np.where(((vis1_2d[i1[3],:]/viscos-1)<=1) == 0))
tbl1100 = np.max(np.where(((vis1_2d[i1[4],:]/viscos-1)<=1) == 0))
tbl1110 = np.max(np.where(((vis1_2d[i1[5],:]/viscos-1)<=1) == 0))
tbl1120 = np.max(np.where(((vis1_2d[i1[6],:]/viscos-1)<=1) == 0))
tbl1130 = np.max(np.where(((vis1_2d[i1[7],:]/viscos-1)<=1) == 0))

blt1060 = y_2d[i1[1],tbl1065]
blt1080 = y_2d[i1[2],tbl1080]
blt1090 = y_2d[i1[3],tbl1090]
blt1100 = y_2d[i1[4],tbl1100]
blt1110 = y_2d[i1[5],tbl1110]
blt1120 = y_2d[i1[6],tbl1120]
blt1130 = y_2d[i1[7],tbl1130]

fig17,axes = plt.subplots(2,4,constrained_layout=True)
axes = axes.reshape(-1)
axes[0].plot(-uv_model_2d[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'b-o',ms=4,label='Modeled')
axes[1].plot(-uv_model_2d[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'b-o',ms=4,label='Modeled')
axes[2].plot(-uv_model_2d[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'b-o',ms=4,label='Modeled')
axes[3].plot(-uv_model_2d[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'b-o',ms=4,label='Modeled')
axes[4].plot(-uv_model_2d[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'b-o',ms=4,label='Modeled')
axes[5].plot(-uv_model_2d[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'b-o',ms=4,label='Modeled')
axes[6].plot(-uv_model_2d[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'b-o',ms=4,label='Modeled')
axes[7].plot(-uv_model_2d[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'b-o',ms=4,label='Modeled')
axes[0].plot(-uv_2d[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'r-o',ms=4,label='Resolved')
axes[1].plot(-uv_2d[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'r-o',ms=4,label='Resolved')
axes[2].plot(-uv_2d[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'r-o',ms=4,label='Resolved')
axes[3].plot(-uv_2d[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'r-o',ms=4,label='Resolved')
axes[4].plot(-uv_2d[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'r-o',ms=4,label='Resolved')
axes[5].plot(-uv_2d[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'r-o',ms=4,label='Resolved')
axes[6].plot(-uv_2d[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'r-o',ms=4,label='Resolved')
axes[7].plot(-uv_2d[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'r-o',ms=4,label='Resolved')
axes[1].plot(-x065_off[:,6],x065_off[:,1]-y_2d[i1[1],0],'k-*',ms=4,label='Experimental')
axes[2].plot(-x080_off[:,6],x080_off[:,1]-y_2d[i1[2],0],'k-*',ms=4,label='Experimental')
axes[3].plot(-x090_off[:,6],x090_off[:,1]-y_2d[i1[3],0],'k-*',ms=4,label='Experimental')
axes[4].plot(-x100_off[:,6],x100_off[:,1]-y_2d[i1[4],0],'k-*',ms=4,label='Experimental')
axes[5].plot(-x110_off[:,6],x110_off[:,1]-y_2d[i1[5],0],'k-*',ms=4,label='Experimental')
axes[6].plot(-x120_off[:,6],x120_off[:,1]-y_2d[i1[6],0],'k-*',ms=4,label='Experimental')
axes[7].plot(-x130_off[:,6],x130_off[:,1]-y_2d[i1[7],0],'k-*',ms=4,label='Experimental')
axes[1].plot(np.linspace(np.min(-uv_2d[i1[1],:]),np.max(-uv_2d[i1[1],:]),10),(y_2d[i1[1],tbl065]-y_2d[i1[1],0])*np.ones(10),'k--',ms=4,label='BL thickness')
axes[2].plot(np.linspace(np.min(-uv_2d[i1[2],:]),np.max(-uv_2d[i1[2],:]),10),(y_2d[i1[2],tbl080]-y_2d[i1[2],0])*np.ones(10),'k--',ms=4,label='BL thickness')
axes[3].plot(np.linspace(np.min(-uv_2d[i1[3],:]),np.max(-uv_2d[i1[3],:]),10),(y_2d[i1[3],tbl090]-y_2d[i1[3],0])*np.ones(10),'k--',ms=4,label='BL thickness')
axes[4].plot(np.linspace(np.min(-uv_2d[i1[4],:]),np.max(-uv_2d[i1[4],:]),10),(y_2d[i1[4],tbl100]-y_2d[i1[4],0])*np.ones(10),'k--',ms=4,label='BL thickness')
axes[5].plot(np.linspace(np.min(-uv_2d[i1[5],:]),np.max(-uv_2d[i1[5],:]),10),(y_2d[i1[5],tbl110]-y_2d[i1[5],0])*np.ones(10),'k--',ms=4,label='BL thickness')
axes[6].plot(np.linspace(np.min(-uv_2d[i1[6],:]),np.max(-uv_2d[i1[6],:]),10),(y_2d[i1[6],tbl120]-y_2d[i1[6],0])*np.ones(10),'k--',ms=4,label='BL thickness')
axes[7].plot(np.linspace(np.min(-uv_2d[i1[7],:]),np.max(-uv_2d[i1[7],:]),10),(y_2d[i1[7],tbl130]-y_2d[i1[7],0])*np.ones(10),'k--',ms=4,label='BL thickness')
for ax,xloc in zip(axes,xx):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(title='At X = '+str(xloc)+' m',title_fontsize=14,fontsize=14,loc='best')
    ax.tick_params(labelsize='medium')
    ax.grid(True, linestyle='-.')
    ax.set_ylim([0,0.25])
    if (ax==axes[0] or ax==axes[4]):
        ax.set_ylabel("$y-y_{wall}$")
    if (ax==axes[4] or ax==axes[5] or ax==axes[6] or ax==axes[7]):
        ax.set_xlabel("$-<\overline{u'v'}>$")
fig17.suptitle("Turbulent Shear Stress -$<\overline{u'v'}>$ (Reference [203])",fontsize=22)

fig17a,axes = plt.subplots(1,1,constrained_layout=True)
a1 = axes.contourf(x_2d,y_2d,-uv_model_2d,cmap='jet',levels=50)
cb1 = fig17a.colorbar(a1,ax=axes)
cb1.ax.tick_params(labelsize=22)
a1.ax.tick_params(labelsize=16)
plt.axis([0.6,1.5,0,1])
axes.set_xlabel('X (m)',fontsize=18)
axes.set_ylabel('Y (m)',fontsize=18)
axes.set_title("Modeled Turbulent Shear Stress -$<\overline{u'v'}>$ (Reference [203])",fontsize=24)

fig17b,axes = plt.subplots(1,1,constrained_layout=True)
a1 = axes.contourf(x_2d,y_2d,-uv_2d,cmap='jet',levels=50)
cb1 = fig17b.colorbar(a1,ax=axes)
cb1.ax.tick_params(labelsize=22)
a1.ax.tick_params(labelsize=16)
plt.axis([0.6,1.5,0,1])
axes.set_xlabel('X (m)',fontsize=18)
axes.set_ylabel('Y (m)',fontsize=18)
axes.set_title("Resolved Turbulent Shear Stress -$<\overline{u'v'}>$ (Reference [203])",fontsize=24)

fig17c,axes = plt.subplots(2,4,constrained_layout=True)
axes = axes.reshape(-1)
axes[0].plot(-uv1_model_2d[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'b-o',ms=4,label='Modeled')
axes[1].plot(-uv1_model_2d[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'b-o',ms=4,label='Modeled')
axes[2].plot(-uv1_model_2d[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'b-o',ms=4,label='Modeled')
axes[3].plot(-uv1_model_2d[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'b-o',ms=4,label='Modeled')
axes[4].plot(-uv1_model_2d[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'b-o',ms=4,label='Modeled')
axes[5].plot(-uv1_model_2d[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'b-o',ms=4,label='Modeled')
axes[6].plot(-uv1_model_2d[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'b-o',ms=4,label='Modeled')
axes[7].plot(-uv1_model_2d[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'b-o',ms=4,label='Modeled')
axes[0].plot(-uv1_2d[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'r-o',ms=4,label='Resolved')
axes[1].plot(-uv1_2d[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'r-o',ms=4,label='Resolved')
axes[2].plot(-uv1_2d[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'r-o',ms=4,label='Resolved')
axes[3].plot(-uv1_2d[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'r-o',ms=4,label='Resolved')
axes[4].plot(-uv1_2d[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'r-o',ms=4,label='Resolved')
axes[5].plot(-uv1_2d[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'r-o',ms=4,label='Resolved')
axes[6].plot(-uv1_2d[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'r-o',ms=4,label='Resolved')
axes[7].plot(-uv1_2d[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'r-o',ms=4,label='Resolved')
axes[1].plot(-x065_off[:,6],x065_off[:,1]-y_2d[i1[1],0],'k-*',ms=4,label='Experimental')
axes[2].plot(-x080_off[:,6],x080_off[:,1]-y_2d[i1[2],0],'k-*',ms=4,label='Experimental')
axes[3].plot(-x090_off[:,6],x090_off[:,1]-y_2d[i1[3],0],'k-*',ms=4,label='Experimental')
axes[4].plot(-x100_off[:,6],x100_off[:,1]-y_2d[i1[4],0],'k-*',ms=4,label='Experimental')
axes[5].plot(-x110_off[:,6],x110_off[:,1]-y_2d[i1[5],0],'k-*',ms=4,label='Experimental')
axes[6].plot(-x120_off[:,6],x120_off[:,1]-y_2d[i1[6],0],'k-*',ms=4,label='Experimental')
axes[7].plot(-x130_off[:,6],x130_off[:,1]-y_2d[i1[7],0],'k-*',ms=4,label='Experimental')
axes[1].plot(np.linspace(np.min(-uv1_2d[i1[1],:]),np.max(-uv1_2d[i1[1],:]),10),(y_2d[i1[1],tbl1065]-y_2d[i1[1],0])*np.ones(10),'k--',ms=4,label='BL thickness')
axes[2].plot(np.linspace(np.min(-uv1_2d[i1[2],:]),np.max(-uv1_2d[i1[2],:]),10),(y_2d[i1[2],tbl1080]-y_2d[i1[2],0])*np.ones(10),'k--',ms=4,label='BL thickness')
axes[3].plot(np.linspace(np.min(-uv1_2d[i1[3],:]),np.max(-uv1_2d[i1[3],:]),10),(y_2d[i1[3],tbl1090]-y_2d[i1[3],0])*np.ones(10),'k--',ms=4,label='BL thickness')
axes[4].plot(np.linspace(np.min(-uv1_2d[i1[4],:]),np.max(-uv1_2d[i1[4],:]),10),(y_2d[i1[4],tbl1100]-y_2d[i1[4],0])*np.ones(10),'k--',ms=4,label='BL thickness')
axes[5].plot(np.linspace(np.min(-uv1_2d[i1[5],:]),np.max(-uv1_2d[i1[5],:]),10),(y_2d[i1[5],tbl1110]-y_2d[i1[5],0])*np.ones(10),'k--',ms=4,label='BL thickness')
axes[6].plot(np.linspace(np.min(-uv1_2d[i1[6],:]),np.max(-uv1_2d[i1[6],:]),10),(y_2d[i1[6],tbl1120]-y_2d[i1[6],0])*np.ones(10),'k--',ms=4,label='BL thickness')
axes[7].plot(np.linspace(np.min(-uv1_2d[i1[7],:]),np.max(-uv1_2d[i1[7],:]),10),(y_2d[i1[7],tbl1130]-y_2d[i1[7],0])*np.ones(10),'k--',ms=4,label='BL thickness')
for ax,xloc in zip(axes,xx):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(title='At X = '+str(xloc)+' m',title_fontsize=14,fontsize=14,loc='best')
    ax.tick_params(labelsize='medium')
    ax.grid(True, linestyle='-.')
    ax.set_ylim([0,0.25])
    if (ax==axes[0] or ax==axes[4]):
        ax.set_ylabel("$y-y_{wall}$")
    if (ax==axes[4] or ax==axes[5] or ax==axes[6] or ax==axes[7]):
        ax.set_xlabel("$-<\overline{u'v'}>$")
fig17c.suptitle("Turbulent Shear Stress -$<\overline{u'v'}>$ (Reference [169])",fontsize=22)

fig17d,axes = plt.subplots(1,1,constrained_layout=True)
a1 = axes.contourf(x_2d,y_2d,-uv1_model_2d,cmap='jet',levels=50)
cb1 = fig17d.colorbar(a1,ax=axes)
cb1.ax.tick_params(labelsize=22)
a1.ax.tick_params(labelsize=16)
plt.axis([0.6,1.5,0,1])
axes.set_xlabel('X (m)',fontsize=18)
axes.set_ylabel('Y (m)',fontsize=18)
axes.set_title("Modeled Turbulent Shear Stress -$<\overline{u'v'}>$ (Reference [169])",fontsize=24)

fig17e,axes = plt.subplots(1,1,constrained_layout=True)
a1 = axes.contourf(x_2d,y_2d,-uv1_2d,cmap='jet',levels=50)
cb1 = fig17e.colorbar(a1,ax=axes)
cb1.ax.tick_params(labelsize=22)
a1.ax.tick_params(labelsize=16)
plt.axis([0.6,1.5,0,1])
axes.set_xlabel('X (m)',fontsize=18)
axes.set_ylabel('Y (m)',fontsize=18)
axes.set_title("Resolved Turbulent Shear Stress -$<\overline{u'v'}>$ (Reference [169])",fontsize=24)

# %% Assignment T3 - The Turbulent Viscosity
visc_ratio = visturb_2d/viscos
visc1_ratio = visturb1_2d/viscos

# Visualize Turbulent viscosity ratio
vmin1 = 1
vmax1 = np.max(np.max(visc_ratio))
v1 = np.linspace(vmin1,vmax1,100)
vmin2 = 1
vmax2 = np.max(np.max(visc1_ratio))
v2 = np.linspace(vmin2,vmax2,100)
fig15,axes = plt.subplots(2,1,constrained_layout=True)
a1 = axes[0].contourf(x_2d,y_2d,visc_ratio,cmap='gist_gray',levels=v1,extend='min')
a2 = axes[1].contourf(x_2d,y_2d,visc1_ratio,cmap='gist_gray',levels=v2,extend='min')
a1.cmap.set_under('black')
a2.cmap.set_under('black')
cb1 = fig15.colorbar(a1,ax=axes[0])
cb2 = fig15.colorbar(a2,ax=axes[1])
cb1.ax.tick_params(labelsize=16)
cb2.ax.tick_params(labelsize=16)
a1.ax.tick_params(labelsize=16)
a2.ax.tick_params(labelsize=16)
axes[0].set_title(r'Turbulent Viscosity Ratio - $\nu_{t}/\nu$ (Paper Data [203])',fontsize=22)
axes[1].set_title(r'Turbulent Viscosity Ratio - $\nu_{t}/\nu$ (Journal Data [169])',fontsize=22)
for ax in axes:
    ax.axis([0.6,1.5,0,1])
    ax.set_xlabel('X (m)',fontsize=18)
    ax.set_ylabel('Y (m)',fontsize=18)

fig19,axes = plt.subplots(2,4,constrained_layout=True)
axes = axes.reshape(-1)
axes[0].plot(visc_ratio[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'b-o',ms=4,label='Paper Data [203]')
axes[1].plot(visc_ratio[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'b-o',ms=4,label='Paper Data [203]')
axes[2].plot(visc_ratio[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'b-o',ms=4,label='Paper Data [203]')
axes[3].plot(visc_ratio[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'b-o',ms=4,label='Paper Data [203]')
axes[4].plot(visc_ratio[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'b-o',ms=4,label='Paper Data [203]')
axes[5].plot(visc_ratio[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'b-o',ms=4,label='Paper Data [203]')
axes[6].plot(visc_ratio[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'b-o',ms=4,label='Paper Data [203]')
axes[7].plot(visc_ratio[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'b-o',ms=4,label='Paper Data [203]')
axes[0].plot(visc1_ratio[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'r-o',ms=4,label='Journal Data [169]')
axes[1].plot(visc1_ratio[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'r-o',ms=4,label='Journal Data [169]')
axes[2].plot(visc1_ratio[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'r-o',ms=4,label='Journal Data [169]')
axes[3].plot(visc1_ratio[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'r-o',ms=4,label='Journal Data [169]')
axes[4].plot(visc1_ratio[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'r-o',ms=4,label='Journal Data [169]')
axes[5].plot(visc1_ratio[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'r-o',ms=4,label='Journal Data [169]')
axes[6].plot(visc1_ratio[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'r-o',ms=4,label='Journal Data [169]')
axes[7].plot(visc1_ratio[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'r-o',ms=4,label='Journal Data [169]')
for ax,xloc in zip(axes,xx):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(title='At X = '+str(xloc)+' m',title_fontsize=14,fontsize=14,loc='best')
    ax.tick_params(labelsize='medium')
    ax.grid(True, linestyle='-.')
    if (ax==axes[0] or ax==axes[4]):
        ax.set_ylabel("$y-y_{wall}$")
    if (ax==axes[4] or ax==axes[5] or ax==axes[6] or ax==axes[7]):
        ax.set_xlabel(r"$\nu_{t}/\nu$")
fig19.suptitle(r'Turbulent Viscosity Ratio - $\nu_{t}/\nu$',fontsize=22)

# term 1,5 (i=1 and j=1)
# term 2,6 (i=1 and j=2)
# term 3,7 (i=2 and j=1)
# term 4,8 (i=2 and j=2)
term1,dummy1 = dphidx_dy(x_2d_new,y_2d_new,visturb_2d*dudx)
dummy2,term2 = dphidx_dy(x_2d_new,y_2d_new,visturb_2d*dudy)
term3,dummy3 = dphidx_dy(x_2d_new,y_2d_new,visturb_2d*dvdx)
dummy4,term4 = dphidx_dy(x_2d_new,y_2d_new,visturb_2d*dvdy)

term5,dummy5 = dphidx_dy(x_2d_new,y_2d_new,uu_2d)
dummy6,term6 = dphidx_dy(x_2d_new,y_2d_new,uv_2d)
term7,dummy7 = dphidx_dy(x_2d_new,y_2d_new,uv_2d)
dummy8,term8 = dphidx_dy(x_2d_new,y_2d_new,vv_2d)

div1_i1 = term1 + term2
div1_i2 = term3 + term4
div2_i1 = -term5 - term6
div2_i2 = -term7 - term8

fig20,axes = plt.subplots(2,4,constrained_layout=True)
axes = axes.reshape(-1)
axes[0].plot(div1_i1[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'b-o',ms=4,label=r'Modeled')
axes[1].plot(div1_i1[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'b-o',ms=4,label=r'Modeled')
axes[2].plot(div1_i1[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'b-o',ms=4,label=r'Modeled')
axes[3].plot(div1_i1[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'b-o',ms=4,label=r'Modeled')
axes[4].plot(div1_i1[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'b-o',ms=4,label=r'Modeled')
axes[5].plot(div1_i1[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'b-o',ms=4,label=r'Modeled')
axes[6].plot(div1_i1[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'b-o',ms=4,label=r'Modeled')
axes[7].plot(div1_i1[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'b-o',ms=4,label=r'Modeled')
axes[0].plot(div2_i1[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'r-o',ms=4,label=r"Resolved")
axes[1].plot(div2_i1[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'r-o',ms=4,label=r"Resolved")
axes[2].plot(div2_i1[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'r-o',ms=4,label=r"Resolved")
axes[3].plot(div2_i1[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'r-o',ms=4,label=r"Resolved")
axes[4].plot(div2_i1[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'r-o',ms=4,label=r"Resolved")
axes[5].plot(div2_i1[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'r-o',ms=4,label=r"Resolved")
axes[6].plot(div2_i1[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'r-o',ms=4,label=r"Resolved")
axes[7].plot(div2_i1[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'r-o',ms=4,label=r"Resolved")
for ax,xloc in zip(axes,xx):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(title='At X = '+str(xloc)+' m',title_fontsize=14,fontsize=14,loc='best')
    ax.tick_params(labelsize='medium')
    ax.grid(True, linestyle='-.')
    ax.set_ylim([0,0.2])
    if (ax==axes[0] or ax==axes[4]):
        ax.set_ylabel("$y-y_{wall}$")
fig20.suptitle(r"Modeled ($\frac{\partial}{\partial x_{j}}(\langle \nu_{t} \frac{\partial \overline{v_{1}}}{\partial x_{j}} \rangle)$) and Resolved ($- \frac{\partial \langle v'_{1}v'_{j} \rangle}{\partial x_{j}}$) Shear Stress",fontsize=22)

fig21,axes = plt.subplots(2,4,constrained_layout=True)
axes = axes.reshape(-1)
axes[0].plot(div1_i2[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'b-o',ms=4,label=r'Modeled')
axes[1].plot(div1_i2[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'b-o',ms=4,label=r'Modeled')
axes[2].plot(div1_i2[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'b-o',ms=4,label=r'Modeled')
axes[3].plot(div1_i2[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'b-o',ms=4,label=r'Modeled')
axes[4].plot(div1_i2[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'b-o',ms=4,label=r'Modeled')
axes[5].plot(div1_i2[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'b-o',ms=4,label=r'Modeled')
axes[6].plot(div1_i2[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'b-o',ms=4,label=r'Modeled')
axes[7].plot(div1_i2[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'b-o',ms=4,label=r'Modeled')
axes[0].plot(div2_i2[i1[0],:],y_2d[i1[0],:]-y_2d[i1[0],0],'r-o',ms=4,label=r"Resolved")
axes[1].plot(div2_i2[i1[1],:],y_2d[i1[1],:]-y_2d[i1[1],0],'r-o',ms=4,label=r"Resolved")
axes[2].plot(div2_i2[i1[2],:],y_2d[i1[2],:]-y_2d[i1[2],0],'r-o',ms=4,label=r"Resolved")
axes[3].plot(div2_i2[i1[3],:],y_2d[i1[3],:]-y_2d[i1[3],0],'r-o',ms=4,label=r"Resolved")
axes[4].plot(div2_i2[i1[4],:],y_2d[i1[4],:]-y_2d[i1[4],0],'r-o',ms=4,label=r"Resolved")
axes[5].plot(div2_i2[i1[5],:],y_2d[i1[5],:]-y_2d[i1[5],0],'r-o',ms=4,label=r"Resolved")
axes[6].plot(div2_i2[i1[6],:],y_2d[i1[6],:]-y_2d[i1[6],0],'r-o',ms=4,label=r"Resolved")
axes[7].plot(div2_i2[i1[7],:],y_2d[i1[7],:]-y_2d[i1[7],0],'r-o',ms=4,label=r"Resolved")
for ax,xloc in zip(axes,xx):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(title='At X = '+str(xloc)+' m',title_fontsize=14,fontsize=14,loc='best')
    ax.tick_params(labelsize='medium')
    ax.grid(True, linestyle='-.')
    ax.set_ylim([0,0.2])
    if (ax==axes[0] or ax==axes[4]):
        ax.set_ylabel("$y-y_{wall}$")
fig21.suptitle(r"Modeled($\frac{\partial}{\partial x_{j}}(\langle \nu_{t} \frac{\partial \overline{v_{2}}}{\partial x_{j}} \rangle)$) and Resolved ($- \frac{\partial \langle v'_{2}v'_{j} \rangle}{\partial x_{j}}$) Shear Stress",fontsize=22)

# %% For visualization - location of experimental points
fig22 = plt.figure()
axes = fig22.add_axes([0.1,0.1,0.8,0.8])
axes.plot(x_2d[:,0],y_2d[:,0],lw=4,color='black')
axes.plot(x_2d[0,:],y_2d[0,:],lw=4,color='black')
axes.plot(x_2d[:,-1],y_2d[:,-1],lw=4,color='black')
axes.plot(x_2d[-1,:],y_2d[-1,:],lw=4,color='black')
axes.plot(x065_off[:,0],np.linspace(x065_off[-1,1],0.6,x065_off.shape[0]),'b--',label='X = 0.65 m')
axes.plot(x080_off[:,0],np.linspace(x080_off[-1,1],0.6,x080_off.shape[0]),'g--',label='X = 0.80 m')
axes.plot(x090_off[:,0],np.linspace(x090_off[-1,1],0.6,x090_off.shape[0]),'r--',label='X = 0.90 m')
axes.plot(x100_off[:,0],np.linspace(x100_off[-1,1],0.6,x100_off.shape[0]),'c--',label='X = 0.10 m')
axes.plot(x110_off[:,0],np.linspace(x110_off[-1,1],0.6,x110_off.shape[0]),'m--',label='X = 0.11 m')
axes.plot(x120_off[:,0],np.linspace(x120_off[-1,1],0.6,x120_off.shape[0]),'y--',label='X = 0.12 m')
axes.plot(x130_off[:,0],np.linspace(x130_off[-1,1],0.6,x130_off.shape[0]),'k--',label='X = 0.13 m')
axes.legend(fontsize=14,loc='best')
axes.tick_params(axis="x", labelsize=15)
axes.tick_params(axis="y", labelsize=15)
axes.set_xlabel('x [m]',fontsize=20)
axes.set_ylabel('y [m]',fontsize=20)
axes.set_title('Location of Points for Experimental Data',fontsize=24)