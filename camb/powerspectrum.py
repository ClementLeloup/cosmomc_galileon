import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

def read(name):
    l = np.genfromtxt(name, usecols=(0))
    Cl = np.genfromtxt(name, usecols=(1))
    return l, Cl
    

liste = ['../../../CosmoMC2015/cosmomc/mycamb/galileon-camb/results/Cls/lcdm/lcdm_scalCovCls.dat', '../../../CosmoMC2015/cosmomc_galileon/camb/test_scalCovCls.dat', 'test_scalCovCls.dat']


l1, Cl1, = read(liste[1])
l2, Cl2, = read(liste[2])
l_lcdm, Cl_lcdm, = read(liste[0])

planck = pyfits.open('../../../CosmoMC2015/cosmomc/mycamb/galileon-camb/COM_PowerSpect_CMB_R2.01.fits')
lowl = planck[1].data
highl = planck[2].data

fig = plt.figure(figsize=(12 ,12))

gs1 = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
#gs1.update(left=0.10, right=0.45, wspace=0.05)
ax1 = plt.subplot(gs1[0])
ax2 = plt.subplot(gs1[1])

maskerror = lowl['ELL']<32
reducedlowl = lowl['ELL'][maskerror]
reducedlowCl = lowl['D_ELL'][maskerror]
reducederrdown = lowl['ERRDOWN'][maskerror]
reducederrup = lowl['ERRUP'][maskerror]

ax1.semilogx()
ax1.plot(l1, Cl1, 'r', label=r'Uncoupled Galileon model best fit Neveu et al. (2016)')
ax1.plot(l2, Cl2, 'g', label=r'Uncoupled Galileon model best fit Neveu et al. (2016)')
ax1.plot(l_lcdm, Cl_lcdm, 'k', linestyle='--', label=r'$\Lambda CDM$')
ax1.errorbar(reducedlowl, reducedlowCl, color='royalblue', ecolor='royalblue', yerr=[reducederrdown, reducederrup], fmt='o')
ax1.errorbar(highl['ELL'], highl['D_ELL'], color='royalblue', ecolor='royalblue', yerr=highl['ERR'], fmt='o', label=r'Planck 2015')
ax1.set_ylim([0, 10000])
ax1.set_xlim([2, 2000])
#ax1.text(80, 700, r'$Galileon \ 2$', fontsize=20)
ax1.legend(loc="upper right", numpoints=1)
ax1.set_ylabel(r'$\ell(\ell+1)\frac{C_{\ell}^{TT}}{2\pi} \ [\mu K^{2}]$', fontsize=17)
ax1.set_xlabel(r'$\ell$', fontsize=17)

ax2.semilogx()
ax2.plot(l1,Cl2/Cl1-1)
#ax2.plot([2, 2000], [0.01, 0.01], linestyle=':', color='k')
ax2.plot([2, 2000], [0., 0.], linestyle='--', color='k')
#ax2.plot([2, 2000], [-0.01, -0.01], linestyle=':', color='k')
ax2.set_ylim([-0.01, 0.01])
ax2.set_xlim([2, 2000])
ax2.set_ylabel(r'$1-\frac{C_{\ell}^{TT}(new\ CAMB)}{C_{\ell}^{TT}(old\ CAMB)}$', fontsize=17)
ax2.set_xlabel(r'$\ell$', fontsize=17)

#fig.savefig("../plot/powerspectrum_jeremy_uncoupled_30000pts.png")
plt.show()
