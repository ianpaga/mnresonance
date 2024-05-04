import matplotlib
import matplotlib.pyplot as plt
import numpy as np 

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['lines.linewidth'] = 1.5

xslices=1
yslices=1
costhslices=3000 

#dtheta=np.pi/300
#dtheta=1.0
dtheta=2.0/costhslices


time = np.genfromtxt("time.dat",delimiter='')

rhoex = np.genfromtxt("rhoex.dat",delimiter='')
rhoexbar = np.genfromtxt("rhoexbar.dat",delimiter='')

rhoee = np.genfromtxt("rhoee.dat",delimiter='')
rhoeebar = np.genfromtxt("rhoeebar.dat",delimiter='')

mnr_nu = np.genfromtxt("mnr_nu.dat",delimiter='')
mnr_nu0 = np.genfromtxt("mnr_nu0.dat",delimiter='')
mnr_lam = np.genfromtxt("mnr_lam.dat",delimiter='')
mnr_vac = np.genfromtxt("mnr_vac.dat",delimiter='')

data_v = np.genfromtxt("test_v0.dat",delimiter="")
data_v_init = np.genfromtxt("test_v0_init.dat",delimiter="")
data_v_nXXX = np.genfromtxt("test_v0_nXXX.dat",delimiter="")


plt.clf()
fig, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,figsize=(12, 11))

ax1.semilogy(time, rhoex*dtheta,label=r'$\rho_{ex}$',color='blue')
ax1.semilogy(time, rhoexbar*dtheta, label=r'$\bar{\rho}_{ex}$',color='red')
ax1.set_xlabel('Time[s]')
ax1.set_ylabel(r'$\rho_{ex},\bar{\rho}_{ex}$')
ax1.legend(loc=1)

ax2.plot(time, rhoee/rhoee[0],label=r'$\rho_{ee}/\rho_{ee}(t=0)$', color='blue')
ax2.plot(time, rhoeebar/rhoeebar[0], label=r'$\bar{\rho}_{ee}/\bar{\rho}_{ee}(t=0)$', color='red')
ax2.set_xlabel('Time[s]')
ax2.set_ylabel(r'$\rho_{ee},\bar{\rho}_{ee}$')
ax2.legend(loc=1)

#ax3.plot(time, mnr_vac, color='purple', label=r'$\omega=%.2f \ \mathrm{km^{-1}}$' %mnr_vac[0])
ax3.plot(time, mnr_lam, color='green', label=r'$\lambda$')
ax3.plot(time, abs(mnr_nu), color='orange', label=r'$|H_{ee}^{\nu\nu}-H_{xx}^{\nu\nu}|_{\cos\theta_0^+}$')
ax3.plot(time, abs(mnr_nu0), color='red', label=r'$|H_{ee}^{\nu\nu}-H_{xx}^{\nu\nu}|_{\cos\theta_0^-}$')
#ax3.set_yscale('symlog') 
ax3.set_xlabel('Time[s]')
ax3.set_ylabel(r'Potentials $\mathrm{[km^{-1}}]$')
ax3.legend(loc=0)

ax4.plot(time, mnr_nu + mnr_lam, color='black', label=r'$(H_{ee}^{\nu\nu}-H_{xx}^{\nu\nu})_{\cos\theta_0^+} + \lambda$')
ax4.plot(time, mnr_nu0 + mnr_lam, color='blue', label=r'$(H_{ee}^{\nu\nu}-H_{xx}^{\nu\nu})_{\cos\theta_0^-} + \lambda$')
#ax4.plot(time, mnr_nu + mnr_lam, color='black', label=r'$\mu + \lambda$')
ax4.axhline(y=0,linestyle='dashed',color='red')
ax4.set_xlabel('Time[s]')
ax4.set_ylabel('MNR condition $\mathrm{[km^{-1}}]$')
ax4.legend(loc=0)

#*#* begin angular distributions #*#*

costhmin= -1.0e0
costhmax= +1.0e0

costharray=np.zeros(costhslices)
cosdtharray=np.zeros(costhslices)

test_rhoee = np.zeros(costhslices);test_rhoxx = np.zeros(costhslices);test_rhoeebar = np.zeros(costhslices);
test_rhoxxbar = np.zeros(costhslices);test_realrhoex = np.zeros(costhslices);test_imagrhoex = np.zeros(costhslices);
test_realrhoexbar = np.zeros(costhslices);test_imagrhoexbar = np.zeros(costhslices);

for k in range(costhslices):
    costharray[k] = costhmin + ((costhmax-costhmin)*(0.5e0+float(k)))*1.0e0/costhslices
    cosdtharray[k] = (costhmax-costhmin)*1.0e0/(costhslices)
    
# initial
for i in range(0, 1):
    for j in range(0,1):
        for k in range(costhslices):
            test_rhoee[k] = data_v_init[i*yslices*costhslices*8+j*costhslices*8+k*8+1]
            test_rhoeebar[k]  = data_v_init[i*yslices*costhslices*8+j*costhslices*8+k*8+5]

ax5.plot( costharray, test_rhoee-test_rhoeebar ,color='green', linestyle='dashed', label=r'$g(v)$')
ax5.axhline(y=0,color='black', linestyle='dashed')
ax5.plot( costharray, test_rhoee ,color='red',linestyle='dashed', label=r'$\nu_e(t=0)$')
ax5.plot( costharray, test_rhoeebar ,color='blue',linestyle='dashed', label=r'$\bar{\nu}_e(t=0)$')

# middle n500
for i in range(0, 1):
    for j in range(0,1):
        for k in range(costhslices):
            test_rhoee[k] = data_v_nXXX[i*yslices*costhslices*8+j*costhslices*8+k*8+1]
            test_rhoeebar[k]  = data_v_nXXX[i*yslices*costhslices*8+j*costhslices*8+k*8+5]

#ax5.plot( costharray, test_rhoee ,color='red',label=r'$\nu_e$')
#ax5.plot( costharray, test_rhoeebar ,color='blue',label=r'$\bar{\nu}_e$')

ax5.set_xlabel(r'$ v = \cos{\theta} $')
ax5.set_ylabel(r'$\rho_{ee}(\cos\theta),\bar{\rho}_{ee}(\cos\theta)$')
ax5.title.set_text(r'$t=t^{\prime}$')
ax5.legend()

# initial
for i in range(0, 1):
    for j in range(0,1):
        for k in range(costhslices):
            test_rhoee[k] = data_v_init[i*yslices*costhslices*8+j*costhslices*8+k*8+1]
            test_rhoeebar[k]  = data_v_init[i*yslices*costhslices*8+j*costhslices*8+k*8+5]

ax6.plot( costharray, test_rhoee ,color='red',linestyle='dashed', label=r'$\nu_e(t=0)$')
ax6.plot( costharray, test_rhoeebar ,color='blue',linestyle='dashed', label=r'$\bar{\nu}_e(t=0)$')

# final
for i in range(0, 1):
    for j in range(0,1):
        for k in range(costhslices):
            test_rhoee[k] = data_v[i*yslices*costhslices*8+j*costhslices*8+k*8+1]
            test_rhoeebar[k]  = data_v[i*yslices*costhslices*8+j*costhslices*8+k*8+5]

ax6.plot( costharray, test_rhoee ,'r', label=r'$\nu_e$')
ax6.plot( costharray, test_rhoeebar ,'b',label=r'$\bar{\nu}_e$')

ax6.set_xlabel(r'$ v = \cos{\theta} $')
ax6.set_ylabel(r'$\rho_{ee}(\cos\theta),\bar{\rho}_{ee}(\cos\theta)$')
ax6.title.set_text(r'$t=t_f$')
ax6.legend()

#*#* end angular distributions #*#*


plt.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.4)


plt.savefig('fig_mnr.pdf')
#plt.show()
