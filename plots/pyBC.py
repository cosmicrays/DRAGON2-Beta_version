#!/bin/bash/python
"""
###############################################################################
# Licence: DRAGON TEAM                                                        #
# Author: Giuseppe Di Bernardo <g.dibernardo@gmail.com>                       # 
# Rev: Nov.7, 2013                                                            #
###############################################################################
"""
import os
import numpy as np
import math
from matplotlib import pyplot as plt
from matplotlib import rc, rcParams 
from scipy.interpolate import interp1d  # useful for the force-field modulation
import pylab
import astropysics
import pyfits
import atpy                             # high-level accessing table data

# activate latex text rendering 
rc('text', usetex=True)
# Change all fonts to 'Computer Modern'
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('xtick', labelsize=18)
rc('ytick', labelsize=18)
rcParams['legend.numpoints'] = 1

print '##########################################################'
print ' Welcome to DRAGONpy. Today, we analyze the B/C data      '
print '##########################################################'

#/*****************************************************************************
# CONSTANTS: physical constants and conversion factors (mostly in cgs units)
#/*****************************************************************************
from astropysics.constants import c, G, mp, me, kb, h, hbar
from astropysics.constants import ergperev, secperday, secperyr, cmperpc
"""i.e., secperyr is a tropical year in sec, cmperpc is a parsec in cm, etc,"""

sec_yr = pow(secperyr,-1)       # to pass from seconds to year
pi = math.pi                    # pi-greek
pc_mt = cmperpc*1.e-2           # parsec in meters
m_to_pc= pow(pc_mt,-1)          # to pass from meters in parsec
m_p = .939                      # proton mass in GeV/c^2
m_e = 0.000510998918            # electron mass in GeV/c^2
b = 1.4e-16                     # energy loss rate in GeV^-1*s^-1
GeV_erg = ergperev*1.e9         # from GeV to erg
c2kB = 6.50966198               # divided by 10^36

#/*****************************************************************************
# FILENAMES
#/*****************************************************************************
ASTROFOLD = os.path.join(os.environ['HOME'],'Astro','Codes')
OUTFOLD   = os.path.join(os.environ['HOME'],'Astro','Codes','DRAGONpy','outpy')
INFOLD    = os.path.join(os.environ['HOME'],'Astro','Codes','DRAGONpy','inpy')
DATAFOLD  = os.path.join(os.environ['HOME'],'Astro','Codes','DRAGONpy','datapy')
NAME      = 'run_2D_v.2_spectrum.fits.gz' # this choice is up to you 
FILE      = INFOLD + '/' + NAME

norm_bkg = 1                 # normalization factor due to the background
phi_pro  = 0.7               # modulation potential for protons 

#/*****************************************************************************
# COMPUTATION
#/*****************************************************************************
hdulist = pyfits.open(FILE) 
hdulist.info() 
prihdr = hdulist[0].header  # the primary header 
n_ext  = len(hdulist)

"""
Looking at the properties of meta-data stored in columns
"""
table_hdu    = hdulist[0]              
table_header = table_hdu.header     # header attribute of TABLE 
zmin   = table_header['zmin']
zmax   = table_header['zmax']
rmin   = table_header['rmin']
rmax   = table_header['rmax']
izsun  = table_header['izsun']
irsun  = table_header['irsun']
dimz   = table_header['dimz']
dimr   = table_header['dimr']
robs   = table_header['robs']
zobs   = table_header['zobs']
emin   = table_header['ekmin']
ek_fac = table_header['ekin_fac']
dimE   = table_header['dimE'] 

""" 
Energy scale in GeV units 
"""
E = np.zeros(dimE)
for i in xrange(0,dimE):
	E[i] = emin*pow(ek_fac,i)

Elog = np.zeros(dimE)
for i in xrange(0,dimE):
	Elog[i] = math.log(E[i])

for i in xrange(1,n_ext):
	if hdulist[i].header['A'] == 10 and hdulist[i].header['Z_'] == 5 :
		table_b10 = hdulist[i].data
	elif hdulist[i].header['A'] == 11 and hdulist[i].header['Z_'] == 5 :
		table_b11 = hdulist[i].data
	elif hdulist[i].header['A'] == 12 and hdulist[i].header['Z_'] == 6 :
		table_c12 = hdulist[i].data
	elif hdulist[i].header['A'] == 13 and hdulist[i].header['Z_'] == 6 :
		table_c13 = hdulist[i].data
	elif hdulist[i].header['A'] == 14 and hdulist[i].header['Z_'] == 6 :
		table_c14 = hdulist[i].data		

print "I am modulating boron 10..."
b10 = table_b10
b10_mod = np.zeros(int(dimE)) 
A = ((E + m_p)**2 - m_p**2) / (((E + m_p) + (5./10.)*phi_pro)**2 - m_p**2)
F = interp1d(Elog,b10,bounds_error=True,kind='cubic')
e_shift = np.zeros(dimE)
for i in xrange(dimE-1):
   	e_shift[i] = math.log(E[i] + (5./10.)*phi_pro)
b10_mod = A * F(e_shift)
print 'modulated protons:',b10_mod

print "I am modulating boron 11..."
b11 = table_b11
b11_mod = np.zeros(int(dimE)) 
A = ((E + m_p)**2 - m_p**2) / ((E + m_p + (5./11.)*phi_pro)**2 - m_p**2)
F = interp1d(Elog,b11,bounds_error=True,kind='cubic')
e_shift = np.zeros(dimE)
for i in xrange(dimE-1):
   	e_shift[i] = math.log(E[i] + (5./11.)*phi_pro)
b11_mod = A * F(e_shift)
print 'modulated protons:',b11_mod

print "I am modulating carbon 12..."
c12 = table_c12
c12_mod = np.zeros(int(dimE)) 
A = ((E + m_p)**2 - m_p**2) / (((E + m_p) + (6./12.)*phi_pro)**2 - m_p**2)
F = interp1d(Elog,c12,bounds_error=True,kind='cubic')
e_shift = np.zeros(dimE)
for i in xrange(dimE-1):
   	e_shift[i] = math.log(E[i] + (6./12.)*phi_pro)
c12_mod = A * F(e_shift)
print 'modulated carbon 12:',c12_mod

print "I am modulating carbon 13..."
c13 = table_c13
c13_mod = np.zeros(int(dimE)) 
A = ((E + m_p)**2 - m_p**2) / ((E + m_p + (6./13.)*phi_pro)**2 - m_p**2)
F = interp1d(Elog,c13,bounds_error=True,kind='cubic')
e_shift = np.zeros(dimE)
for i in xrange(dimE-1):
   	e_shift[i] = math.log(E[i] + (6./13.)*phi_pro)
c13_mod = A * F(e_shift)
print 'modulated carbon 13:',c13_mod

print "I am modulating carbon 14..."
c14 = table_c14
c14_mod = np.zeros(int(dimE)) 
A = ((E + m_p)**2 - m_p**2) / ((E + m_p + (6./14.)*phi_pro)**2 - m_p**2)
F = interp1d(Elog,c14,bounds_error=False)
e_shift = np.zeros(dimE)
for i in xrange(dimE-1):
   	e_shift[i] = math.log(E[i] + (6./14.)*phi_pro)
c14_mod = A * F(e_shift)
print 'modulated carbon 14:',c14_mod

boron = (table_b10 + table_b11)
carbon = (table_c12 + table_c13 + table_c14)
borontocarbon = boron/carbon 
B_spec   = (b10_mod + b11_mod)           # modulated boron spectrum 
C_spec   = (c12_mod + c13_mod + c14_mod) # modulated carbon spectrum 
B_C      = (B_spec/C_spec)               # modulated boron to carbon ratio

#/*****************************************************************************
# HANDLING EXPERIMENTAL DATA
#/*****************************************************************************
data = DATAFOLD + '/' + 'HEAO3_BC.txt'
Emin,Emax = np.loadtxt(data,skiprows=2,usecols=(1,2),unpack=True)
Emean     = np.loadtxt(data,skiprows=2,usecols=(0,),unpack=True)
Ratio     = np.loadtxt(data,skiprows=2,usecols=(3,),unpack=True)
Ratio_low = np.loadtxt(data,skiprows=2,usecols=(8,),unpack=True)
Ratio_up  = np.loadtxt(data,skiprows=2,usecols=(9,),unpack=True)

#/*****************************************************************************
# CURVE PLOTTING
#/*****************************************************************************
print 'Please wait...I am plotting for you'
plt.plot(E, borontocarbon, ls='--' , color='black',lw=2)
plt.plot(E, B_C, ls='-' , color='blue',lw=2)
pam, = plt.plot(Emean,Ratio,'ro')
plt.errorbar(Emean,Ratio,xerr=[Emin, Emax],
   	         yerr=[Ratio_low,Ratio_up],fmt=None, ecolor='b')
 
plt.axis([1.e-1,5.e2, 0.,0.4],interpolation='none')
plt.xlabel(r'E [GeV/nuc]',fontsize=20)
plt.ylabel(r'B/C',fontsize=20)
plt.xscale('log')
#plt.yscale('log')
plt.legend([pam],[r'\textbf{\textit{HEAO3}}'],loc='upper right')
plt.savefig(OUTFOLD + '/' + 'pyBC.eps',format='eps',bbox_inches='tight', dpi=100)
plt.show()

exit ()
