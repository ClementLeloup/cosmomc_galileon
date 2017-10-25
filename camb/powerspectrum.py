import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

def read(name):
    l = np.genfromtxt(name, usecols=(0))
    Cl = np.genfromtxt(name, usecols=(1))
    return l, Cl
    

liste = ['../../cosmomc/mycamb/galileon-camb/results/Cls/lcdm/lcdm_scalCovCls.dat', 'test_scalCovCls.dat']


l1, Cl1, = read(liste[1])
l_lcdm, Cl_lcdm, = read(liste[0])

planck = pyfits.open('../../cosmomc/mycamb/galileon-camb/COM_PowerSpect_CMB_R2.01.fits')
lowl = planck[1].data
highl = planck[2].data

fig = plt.figure(figsize=(12 ,9))

maskerror = lowl['ELL']<32
reducedlowl = lowl['ELL'][maskerror]
reducedlowCl = lowl['D_ELL'][maskerror]
reducederrdown = lowl['ERRDOWN'][maskerror]
reducederrup = lowl['ERRUP'][maskerror]

plt.semilogx()
plt.plot(l1, Cl1, 'r', label=r'Uncoupled Galileon model best fit Neveu et al. (2016)')
plt.plot(l_lcdm, Cl_lcdm, 'k', linestyle='--', label=r'$\Lambda CDM$')
plt.errorbar(reducedlowl, reducedlowCl, color='royalblue', ecolor='royalblue', yerr=[reducederrdown, reducederrup], fmt='o')
plt.errorbar(highl['ELL'], highl['D_ELL'], color='royalblue', ecolor='royalblue', yerr=highl['ERR'], fmt='o', label=r'Planck 2015')
plt.ylim([0, 10000])
plt.xlim([2, 2000])
#plt.text(80, 700, r'$Galileon \ 2$', fontsize=20)
plt.legend(loc="upper right", numpoints=1)
plt.ylabel(r'$\ell(\ell+1)\frac{C_{\ell}^{TT}}{2\pi} \ [\mu K^{2}]$', fontsize=17)
plt.xlabel(r'$\ell$', fontsize=17)

#fig.savefig("../plot/powerspectrum_jeremy_uncoupled_30000pts.png")
plt.show()
