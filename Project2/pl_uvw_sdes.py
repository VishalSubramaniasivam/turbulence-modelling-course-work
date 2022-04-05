import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
from IPython import display
plt.rcParams.update({'font.size': 22})

# makes sure figures are updated when using ipython
display.clear_output()

dx=0.05
dz=0.025

ni=34
nj=49
nk=34

viscos=1./8000.

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  read v_1 & transform v_1 to a 3D array (file 1)
# V1 Velocity
uvw = sio.loadmat('u1_sdes.mat')
ttu=uvw['u1_sdes']
u3d1= np.reshape(ttu,(nk,nj,ni))

uvw = sio.loadmat('u2_sdes.mat')
ttu=uvw['u2_sdes']
u3d2= np.reshape(ttu,(nk,nj,ni))

uvw = sio.loadmat('u3_sdes.mat')
ttu=uvw['u3_sdes']
u3d3= np.reshape(ttu,(nk,nj,ni))

uvw = sio.loadmat('u4_sdes.mat')
ttu=uvw['u4_sdes']
u3d4= np.reshape(ttu,(nk,nj,ni))

uvw = sio.loadmat('u5_sdes.mat')
ttu=uvw['u5_sdes']
u3d5= np.reshape(ttu,(nk,nj,ni))

uvw = sio.loadmat('u6_sdes.mat')
ttu=uvw['u6_sdes']
u3d6= np.reshape(ttu,(nk,nj,ni))

uvw = sio.loadmat('u7_sdes.mat')
ttv=uvw['u7_sdes']
u3d7= np.reshape(ttv,(nk,nj,ni))

uvw = sio.loadmat('u8_sdes.mat')
ttv=uvw['u8_sdes']
u3d8= np.reshape(ttv,(nk,nj,ni))

# V2 Velocity
uvw = sio.loadmat('v1_sdes.mat')
ttv=uvw['v1_sdes']
v3d1= np.reshape(ttv,(nk,nj,ni))

uvw = sio.loadmat('v2_sdes.mat')
ttv=uvw['v2_sdes']
v3d2= np.reshape(ttv,(nk,nj,ni))

uvw = sio.loadmat('v3_sdes.mat')
ttv=uvw['v3_sdes']
v3d3= np.reshape(ttv,(nk,nj,ni))

uvw = sio.loadmat('v4_sdes.mat')
ttv=uvw['v4_sdes']
v3d4= np.reshape(ttv,(nk,nj,ni))

uvw = sio.loadmat('v5_sdes.mat')
ttv=uvw['v5_sdes']
v3d5= np.reshape(ttv,(nk,nj,ni))

uvw = sio.loadmat('v6_sdes.mat')
ttv=uvw['v6_sdes']
v3d6= np.reshape(ttv,(nk,nj,ni))

uvw = sio.loadmat('v7_sdes.mat')
ttv=uvw['v7_sdes']
v3d7= np.reshape(ttv,(nk,nj,ni))

uvw = sio.loadmat('v8_sdes.mat')
ttv=uvw['v8_sdes']
v3d8= np.reshape(ttv,(nk,nj,ni))

# V3 Velocity
uvw = sio.loadmat('w1_sdes.mat')
ttw=uvw['w1_sdes']
w3d1= np.reshape(ttw,(nk,nj,ni))

uvw = sio.loadmat('w2_sdes.mat')
ttw=uvw['w2_sdes']
w3d2= np.reshape(ttw,(nk,nj,ni))

uvw = sio.loadmat('w3_sdes.mat')
ttw=uvw['w3_sdes']
w3d3= np.reshape(ttw,(nk,nj,ni))

uvw = sio.loadmat('w4_sdes.mat')
ttw=uvw['w4_sdes']
w3d4= np.reshape(ttw,(nk,nj,ni))

uvw = sio.loadmat('w5_sdes.mat')
ttw=uvw['w5_sdes']
w3d5= np.reshape(ttw,(nk,nj,ni))

uvw = sio.loadmat('w6_sdes.mat')
ttw=uvw['w6_sdes']
w3d6= np.reshape(ttw,(nk,nj,ni))

uvw = sio.loadmat('w7_sdes.mat')
ttw=uvw['w7_sdes']
w3d7= np.reshape(ttw,(nk,nj,ni))

uvw = sio.loadmat('w8_sdes.mat')
ttw=uvw['w8_sdes']
w3d8= np.reshape(ttw,(nk,nj,ni))

# Turbulent Kinetic Energy
uvw = sio.loadmat('te1_sdes.mat')
ttke=uvw['te1_sdes']
te3d1= np.reshape(ttke,(nk,nj,ni))

uvw = sio.loadmat('te2_sdes.mat')
ttke=uvw['te2_sdes']
te3d2= np.reshape(ttke,(nk,nj,ni))

uvw = sio.loadmat('te3_sdes.mat')
ttke=uvw['te3_sdes']
te3d3= np.reshape(ttke,(nk,nj,ni))

uvw = sio.loadmat('te4_sdes.mat')
ttke=uvw['te4_sdes']
te3d4= np.reshape(ttke,(nk,nj,ni))

uvw = sio.loadmat('te5_sdes.mat')
ttke=uvw['te5_sdes']
te3d5= np.reshape(ttke,(nk,nj,ni))

uvw = sio.loadmat('te6_sdes.mat')
ttke=uvw['te6_sdes']
te3d6= np.reshape(ttke,(nk,nj,ni))

uvw = sio.loadmat('te7_sdes.mat')
ttke=uvw['te7_sdes']
te3d7= np.reshape(ttke,(nk,nj,ni))

uvw = sio.loadmat('te8_sdes.mat')
ttke=uvw['te8_sdes']
te3d8= np.reshape(ttke,(nk,nj,ni))

# Specific Dissipation
uvw = sio.loadmat('om1_sdes.mat')
ttom=uvw['om1_sdes']
om3d1= np.reshape(ttom,(nk,nj,ni))

uvw = sio.loadmat('om2_sdes.mat')
ttom=uvw['om2_sdes']
om3d2= np.reshape(ttom,(nk,nj,ni))

uvw = sio.loadmat('om3_sdes.mat')
ttom=uvw['om3_sdes']
om3d3= np.reshape(ttom,(nk,nj,ni))

uvw = sio.loadmat('om4_sdes.mat')
ttom=uvw['om4_sdes']
om3d4= np.reshape(ttom,(nk,nj,ni))

uvw = sio.loadmat('om5_sdes.mat')
ttom=uvw['om5_sdes']
om3d5= np.reshape(ttom,(nk,nj,ni))

uvw = sio.loadmat('om6_sdes.mat')
ttom=uvw['om6_sdes']
om3d6= np.reshape(ttom,(nk,nj,ni))

uvw = sio.loadmat('om7_sdes.mat')
ttom=uvw['om7_sdes']
om3d7= np.reshape(ttom,(nk,nj,ni))

uvw = sio.loadmat('om8_sdes.mat')
ttom=uvw['om8_sdes']
om3d8= np.reshape(ttom,(nk,nj,ni))

# merge 2 files. This means that new ni = 2*ni
u3d=np.concatenate((u3d1,u3d2,u3d3,u3d4,u3d5,u3d6,u3d7,u3d8), axis=0)
v3d=np.concatenate((v3d1,v3d2,v3d3,v3d4,v3d5,v3d6,v3d7,v3d8), axis=0)
w3d=np.concatenate((w3d1,w3d2,w3d3,w3d4,w3d5,w3d6,w3d7,w3d8), axis=0)
te3d=np.concatenate((te3d1,te3d2,te3d3,te3d4,te3d5,te3d6,te3d7,te3d8), axis=0)
om3d=np.concatenate((om3d1,om3d2,om3d3,om3d4,om3d5,om3d6,om3d7,om3d8), axis=0)

# x coordinate direction = index 0, first index
# y coordinate direction = index 1, second index
# z coordinate direction = index 2, third index

ni=len(u3d)

x=dx*ni
z=dz*nk

umean=np.mean(u3d, axis=(0,2))
vmean=np.mean(v3d, axis=(0,2))
wmean=np.mean(w3d, axis=(0,2))
temean=np.mean(te3d, axis=(0,2))
ommean=np.mean(om3d, axis=(0,2))

# face coordinates
yc = np.loadtxt("yc.dat")

# cell cener coordinates
y= np.zeros(nj)
dy=np.zeros(nj)
for j in range (1,nj-1):
# dy = cell width
   dy[j]=yc[j]-yc[j-1]
   y[j]=0.5*(yc[j]+yc[j-1])

y[nj-1]=yc[nj-1]
tauw=viscos*umean[1]/y[1]
ustar=tauw**0.5
yplus=y*ustar/viscos

DNS=np.genfromtxt("LM_Channel_5200_mean_prof.dat", dtype=None,comments="%")
y_DNS=DNS[:,0]
yplus_DNS=DNS[:,1]
u_DNS=DNS[:,2]

DNS=np.genfromtxt("LM_Channel_5200_vel_fluc_prof.dat", dtype=None,comments="%")
u2_DNS=DNS[:,2]
v2_DNS=DNS[:,3]
w2_DNS=DNS[:,4]
uv_DNS=DNS[:,5]

k_DNS=0.5*(u2_DNS+v2_DNS+w2_DNS)

# find equi.distant DNS cells in log-scale
xx=0.
jDNS=[1]*40
for i in range (0,40):
   i1 = (np.abs(10.**xx-yplus_DNS)).argmin()
   jDNS[i]=int(i1)
   xx=xx+0.2

# %% Assignment S2 - Mean Velocity Profile
fig1,axes = plt.subplots(1,1,constrained_layout=True)
axes.semilogx(y,umean/ustar,'b--',label='DES Data')
axes.semilogx(y_DNS[jDNS],u_DNS[jDNS],'r-o',label='DNS Data')
axes.semilogx(0.0175*np.ones(20),np.linspace(0,25,20),'k--',label=r'$y/\delta$ = 0.0175')
axes.set_xlabel(r"$y/\delta$")
axes.set_ylabel(r"$U^+$")
axes.legend(fontsize=16,loc='best')
axes.tick_params(labelsize='medium')
axes.grid(True, linestyle='-.')
fig1.suptitle(r'Time Averaged Velocity ($<\overline{{v_{1}}^{\prime}}>$)',fontsize=22)

# %% Assignment S3 - Resolved Turbulent Shear Stress
uv_mean = np.mean((u3d-umean[None,:,None])*(v3d-vmean[None,:,None]),axis=(0,2))

fig2,axes = plt.subplots(1,1,constrained_layout=True)
axes.semilogx(y,uv_mean,'b--',label='DES Data - Resolved Turbulent Shear Stress')
axes.semilogx(y_DNS[jDNS],uv_DNS[jDNS],'r-o',label='DNS Data - Total Turbulent Shear Stress')
axes.semilogx(0.0175*np.ones(20),np.linspace(-1,0,20),'k--',label=r'$y/\delta$ = 0.0175')
axes.set_xlabel(r"$y/\delta$")
axes.set_ylabel(r'Shear Stress $<{\overline{{v_{1}}^{\prime}}}{\overline{{v_{2}}^{\prime}}}>$')
axes.legend(fontsize=16,loc='best')
axes.tick_params(labelsize='medium')
axes.grid(True, linestyle='-.')
fig2.suptitle(r'Resolved Turbulent Shear Stress',fontsize=22)

# %% Assignment S4 - Turbulent Kinetic Energy
uu_mean = np.mean((u3d-umean[None,:,None])**2,axis=(0,2))
vv_mean = np.mean((v3d-vmean[None,:,None])**2,axis=(0,2))
ww_mean = np.mean((w3d-wmean[None,:,None])**2,axis=(0,2))

#uu_mean = np.mean(u3d*u3d,axis=(0,2)) - (umean*umean)
#vv_mean = np.mean(v3d*v3d,axis=(0,2)) - (vmean*vmean)
#ww_mean = np.mean(w3d*w3d,axis=(0,2)) - (wmean*wmean)

ke_res = 0.5*(uu_mean+vv_mean+ww_mean)
ke_res_DNS =  0.5*(u2_DNS+v2_DNS+w2_DNS)

fig3,axes = plt.subplots(1,1,constrained_layout=True)
axes.semilogx(y,temean,'b--',label='DES Data - Modeled Turbulent Kinetic Energy')
axes.semilogx(y,ke_res,'r-o',label='DES Data - Resolved Turbulent Kinetic Energy')
axes.semilogx(y,ke_res+temean,'k-o',label='DES Data - Total Turbulent Kinetic Energy')
axes.semilogx(y_DNS,ke_res_DNS,'g-o',label='DNS Data - Turbulent Kinetic Energy')
axes.semilogx(0.0175*np.ones(10),np.linspace(0,6,10),'k--',label=r'$y/\delta$ = 0.0175')
axes.set_xlabel(r"$y/\delta$")
axes.set_ylabel(r'TKE')
axes.legend(fontsize=16,loc='best')
axes.tick_params(labelsize='medium')
axes.grid(True, linestyle='-.')
fig3.suptitle(r'Turbulent Kinetic Energy',fontsize=22)

# %% Assignment S5 - Modeled Turbulent Shear Stress
turb_visc = te3d/om3d
dudx,dudy,dudz = np.gradient(u3d,dx,y,dz)
dvdx,dvdy,dvdz = np.gradient(v3d,dx,y,dz)
dwdx,dwdy,dwdz = np.gradient(w3d,dx,y,dz)
tau_12 = -turb_visc * (dudy+dvdx)
tau_12_mean = np.mean(tau_12,axis=(0,2))
tau_12_mean[0] = 0

fig4,axes = plt.subplots(1,1,constrained_layout=True)
axes.semilogx(y,uv_mean,'b--',label='DES Data - Resolved Turbulent Shear Stress')
axes.semilogx(y,tau_12_mean,'g--',label='DES Data - Modeled Turbulent Shear Stress')
axes.semilogx(y_DNS[jDNS],uv_DNS[jDNS],'r-o',label='DNS Data - Total Turbulent Shear Stress')
axes.semilogx(0.0175*np.ones(20),np.linspace(-1,0,20),'k--',label=r'$y/\delta$ = 0.0175')
axes.set_xlabel(r"$y/\delta$")
axes.set_ylabel(r'Shear Stress $<{\overline{{v_{1}}^{\prime}}}{\overline{{v_{2}}^{\prime}}}>$')
axes.legend(fontsize=16,loc='best')
axes.tick_params(labelsize='medium')
axes.grid(True, linestyle='-.')
fig4.suptitle(r'Comparision of Shear Stress',fontsize=22)

fig5,axes = plt.subplots(1,1,constrained_layout=True)
axes.semilogx(y,uv_mean+tau_12_mean,'b-o',label='DES Data')
axes.semilogx(y_DNS[jDNS],uv_DNS[jDNS],'r-o',label='DNS Data')
axes.semilogx(0.0175*np.ones(20),np.linspace(-1,0,20),'k--',label=r'$y/\delta$ = 0.0175')
axes.set_xlabel(r"$y/\delta$")
axes.set_ylabel(r'Shear Stress $<{\overline{{v_{1}}^{\prime}}}{\overline{{v_{2}}^{\prime}}}>$')
axes.legend(fontsize=16,loc='best')
axes.tick_params(labelsize='medium')
axes.grid(True, linestyle='-.')
fig5.suptitle(r'Total Turbulent Shear Stress',fontsize=22)

# %% Assignment S6 - Location of Interface in DES and DDES
C_DES_SA = 0.65
delta_SA = C_DES_SA * np.maximum(dx*np.ones([49,1]),dy.reshape(49,1),dz*np.ones([49,1]))

Int_SA = y[(y.reshape(49,1)>=delta_SA).argmax()]
print('Switch for SA-DES model occurs at y+ = '+str(np.around(Int_SA,3)))

C_DES_SST = 0.61
beta_star = 0.09
delta_SST = C_DES_SST * np.maximum(dx*np.ones([49,1]),dy.reshape(49,1),dz*np.ones([49,1]))
Lt = ((temean**0.5)/ommean)/beta_star
term1 = Lt.reshape(49,1)/(delta_SST)
term1[0] = 1e-10

Int_SST = y[(np.ones([49,1])>=term1).argmin()]
print('Switch for SST-DES model occurs at y+ = '+str(np.around(Int_SST,3)))

term2 = 2*(Lt[1:]/y[1:])
term3 = (500*viscos)/(ommean[1:]*(y[1:]**2))
eta = np.insert(np.maximum(term2,term3),0,[1e-10])
F2 = np.tanh(eta**2)
F_DDES = np.maximum(term1*(1-F2.reshape(49,1)),1)
Int_DDES = y[(np.ones([49,1])>=(term1*(1-F2.reshape(49,1)))).argmin()]
print('Switch for DDES model occurs at y+ = '+str(np.around(Int_DDES,3)))

fig6,axes = plt.subplots(1,1,constrained_layout=True)
axes.plot(umean,y,label=r'DES data : Mean velocity profile $<\overline{u}>$')
axes.plot(u_DNS,y_DNS,label=r'DNS data : Mean velocity profile $<\overline{u}>$')
axes.plot(np.linspace(0,30,20),0.0175*np.ones(20),'k--',lw=2,label=r'S-DES model')
axes.plot(np.linspace(0,30,20),Int_SA*np.ones(20),'b--',lw=2,label=r'SA-DES model')
axes.plot(np.linspace(0,30,20),Int_SST*np.ones(20),'y--',lw=2,label=r'SST-DES model')
axes.plot(np.linspace(0,30,20),Int_DDES*np.ones(20),'r--',lw=2,label=r'DDES model')
axes.set_ylabel(r"$y/\delta$")
axes.legend(fontsize=16,loc='best')
axes.tick_params(labelsize='medium')
axes.grid(True, linestyle='-.')
fig6.suptitle(r'Location of Interfaces for various models',fontsize=22)

# %% Assignment S7 - SAS Turbulent Length Scale
# 1D Length Scale
dua_dy = np.gradient(umean,y)
dua2_dy = np.gradient(dua_dy,y)
DNSdua_dy = np.gradient(u_DNS,y_DNS)
DNSdua2_dy = np.gradient(DNSdua_dy,y_DNS)
# Von-karman constant
kappa = 0.4

LVK_1D = kappa * np.abs(dua_dy/dua2_dy)
DNSLVK_1D = kappa * np.abs(DNSdua_dy/DNSdua2_dy)

# 3D Length Scale
du2_dxdx,du2_dxdy,du2_dxdz = np.gradient(dudx,dx,y,dz)
dv2_dxdx,dv2_dxdy,dv2_dxdz = np.gradient(dvdx,dx,y,dz)
dw2_dxdx,dw2_dxdy,dw2_dxdz = np.gradient(dwdx,dx,y,dz)

du2_dydx,du2_dydy,du2_dydz = np.gradient(dudy,dx,y,dz)
dv2_dydx,dv2_dydy,dv2_dydz = np.gradient(dvdy,dx,y,dz)
dw2_dydx,dw2_dydy,dw2_dydz = np.gradient(dwdy,dx,y,dz)

du2_dzdx,du2_dzdy,du2_dzdz = np.gradient(dudz,dx,y,dz)
dv2_dzdx,dv2_dzdy,dv2_dzdz = np.gradient(dvdz,dx,y,dz)
dw2_dzdx,dw2_dzdy,dw2_dzdz = np.gradient(dwdz,dx,y,dz)

S11 = 0.5*(2*dudx)
S12 = 0.5*(dudy+dvdx)
S13 = 0.5*(dudz+dwdx)
S22 = 0.5*(2*dvdy)
S23 = 0.5*(dvdz+dwdy)
S33 = 0.5*(2*dwdz)

S = np.sqrt(2*((S11**2)+(S22**2)+(S33**2)+(2*(S12**2))+(2*(S13**2))+(2*(S23**2))))

U1 = np.sqrt((du2_dxdx*du2_dxdx)+(du2_dydy*du2_dydy)+(du2_dzdz*du2_dzdz)+
             (dv2_dxdx*dv2_dxdx)+(dv2_dydy*dv2_dydy)+(dv2_dzdz*dv2_dzdz)+
             (dw2_dxdx*dw2_dxdx)+(dw2_dydy*dw2_dydy)+(dw2_dzdz*dw2_dzdz))
U2 = np.sqrt((du2_dxdx**2)+(du2_dxdy**2)+(du2_dxdz**2)+(du2_dydx**2)+(du2_dydy**2)+
             (du2_dydz**2)+(du2_dzdx**2)+(du2_dzdy**2)+(du2_dzdz**2)+
             (dv2_dxdx**2)+(dv2_dxdy**2)+(dv2_dxdz**2)+(dv2_dydx**2)+(dv2_dydy**2)+
             (dv2_dydz**2)+(dv2_dzdx**2)+(dv2_dzdy**2)+(dv2_dzdz**2)+
             (dw2_dxdx**2)+(dw2_dxdy**2)+(dw2_dxdz**2)+(dw2_dydx**2)+(dw2_dydy**2)+
             (dw2_dydz**2)+(dw2_dzdx**2)+(dw2_dzdy**2)+(dw2_dzdz**2))

LVK_3D = kappa*np.abs(S/U1)
LVK_3D_2 = kappa*np.abs(S/U2)
delta = (dx*dy*dx)**(1/3)

fig7,axes = plt.subplots(1,1,constrained_layout=True)
axes.plot(y[:],LVK_1D[:],'k--',label=r'Von-Karman Length Scale 1D')
axes.plot(y[:],np.mean(LVK_3D[:],axis=(0,2))[:],'r--',label=r'Von-Karman Length Scale 3D')
axes.plot(y[:],np.mean(LVK_3D_2[:],axis=(0,2))[:],'g--',label=r'Von-Karman Length Scale 3D - Method 2')
#axes.plot(y[:],delta[:])
#axes.plot(y[:],dy[:])
axes.set_ylabel(r"Length Scale")
axes.set_xlabel(r"$y/\delta$")
axes.legend(fontsize=16,loc='best')
axes.tick_params(labelsize='medium')
axes.grid(True, linestyle='-.')
fig7.suptitle(r'Von-Karman Length Scale',fontsize=22)

fig7a = plt.figure()
axes1 = fig7a.add_axes([0.1, 0.1, 0.8, 0.8]) # main axes
axes2 = fig7a.add_axes([0.2, 0.55, 0.3, 0.3]) # inset axes
# Larger Figure Axes 1
axes1.plot(y[:],LVK_1D[:],'k-o',label=r'DES data - Von-Karman Length Scale 1D')
axes1.plot(y_DNS[4:],DNSLVK_1D[4:],'r-',label=r'DNS data - Von-Karman Length Scale 1D')
axes1.set_ylabel(r"Length Scale")
axes1.set_xlabel(r"$y/\delta$")
axes1.legend(fontsize=16,loc='upper right')
axes1.tick_params(labelsize='medium')
axes1.grid(True, linestyle='-.')
axes1.set_title(r'Von-Karman Length Scale',fontsize=22)
# Insert Figure Axes 2
axes2.plot(y[:],LVK_1D[:],'k-o',label=r'DES data - Von-Karman Length Scale 1D')
axes2.plot(y_DNS[4:],DNSLVK_1D[4:],'r-',label=r'DNS data - Von-Karman Length Scale 1D')
axes2.set_ylabel(r"Length Scale",fontsize=12)
axes2.set_xlabel(r"$y/\delta$",fontsize=12)
axes2.tick_params(labelsize='x-small')
axes2.grid(True, linestyle='-.')
axes2.set_ylim([-0.01,0.05])
axes2.set_xlim([-0.1,0.2])
axes2.set_title(r'Zoomed near the wall',fontsize=12)

# %% Assignment S8 - Estimate of Resolution
ustar_x = -((dx**2)/4) * (du2_dxdx + dv2_dxdx + dw2_dxdx)
vstar_x = -((dy[None,:,None]**2)/4) * (du2_dxdx + dv2_dxdx + dw2_dxdx)
wstar_x = -((dz**2)/4) * (du2_dxdx + dv2_dxdx + dw2_dxdx)

ustar_y = -((dx**2)/4) * (du2_dydy + dv2_dydy + dw2_dydy)
vstar_y = -((dy[None,:,None]**2)/4) * (du2_dydy + dv2_dydy + dw2_dydy)
wstar_y = -((dz**2)/4) * (du2_dydy + dv2_dydy + dw2_dydy)

ustar_z = -((dx**2)/4) * (du2_dzdz + dv2_dzdz + dw2_dzdz)
vstar_z = -((dy[None,:,None]**2)/4) * (du2_dzdz + dv2_dzdz + dw2_dzdz)
wstar_z = -((dz**2)/4) * (du2_dzdz + dv2_dzdz + dw2_dzdz)

Error_x = (dx*dy*dz) * np.sqrt((np.mean(vstar_x**2,axis=(0,2)))+(np.mean(ustar_x**2,axis=(0,2)))+(np.mean(wstar_x**2,axis=(0,2))))
Error_y = (dx*dy*dz) * np.sqrt((np.mean(vstar_y**2,axis=(0,2)))+(np.mean(ustar_y**2,axis=(0,2)))+(np.mean(wstar_y**2,axis=(0,2))))
Error_z = (dx*dy*dz) * np.sqrt((np.mean(vstar_z**2,axis=(0,2)))+(np.mean(ustar_z**2,axis=(0,2)))+(np.mean(wstar_z**2,axis=(0,2))))

fig8,axes = plt.subplots(1,1,constrained_layout=True)
axes.semilogx(y,Error_x,'r-o',label=r'Anisotropic Error in X-direction')
axes.semilogx(y,Error_y,'k-s',label=r'Anisotropic Error in Y-direction')
axes.semilogx(y,Error_z,'g-*',label=r'Anisotropic Error in Z-direction')
axes.set_xlabel(r"$y/\delta$")
axes.set_ylabel(r"Error")
axes.legend(fontsize=16,loc='best')
axes.tick_params(labelsize='medium')
axes.grid(True, linestyle='-.')
fig8.suptitle(r'Anisotropic Error along Cordinate Directions',fontsize=22)