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
print ' Welcome to DRAGONpy. Today, we analyze the protons data  '
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
NAME      = 'run_2D_spectrum.fits.gz' # this choice is up to you 
FILE      = INFOLD + '/' + NAME
#os.mkdir(outfolder) if you want to make a subfolder under your home folder
norm_bkg = 1.0                 # normalization factor due to the background
norm_ext = 1.0                 # normalization factor 

phi_lep  = 0.30                # modulation potential for leptons 
phi_pro  = 0.35                #   """         """    for protons 


#/*****************************************************************************
# COMPUTATION
#/*****************************************************************************
hdulist = pyfits.open(FILE) 
hdulist.info() 
prihdr = hdulist[0].header  # the primary header 
n_ext  = len(hdulist)
#print prihdr

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
	if hdulist[i].header['A'] == 1 and hdulist[i].header['Z_'] == 1:
		table_pro = hdulist[i].data

print "I am modulating protons..."
pro = table_pro
pro_mod = np.zeros(int(dimE)) 
A = ((E + m_p)**2 - m_p**2) / (((E + m_p) + phi_pro)**2 - m_p**2)
F = interp1d(Elog,pro,bounds_error=True,kind='cubic')
e_shift = np.zeros(dimE)
for i in xrange(dimE-1):
  	e_shift[i] = math.log(E[i] + phi_pro)
pro_mod = A * F(e_shift)
print 'modulated protons:',pro_mod 

pro_spec = norm_bkg * pro_mod     # protons spectrum due to the bkg

#/*****************************************************************************
# HANDLING EXPERIMENTAL DATA
#/*****************************************************************************
data = DATAFOLD + '/' + 'PAMELA_protons.txt'
Emin_pam,Emax_pam = np.loadtxt(data,skiprows=2,usecols=(1,2),unpack=True)
Emean_pam = np.loadtxt(data,skiprows=2,usecols=(0,),unpack=True)
Flux_pam  = np.loadtxt(data,skiprows=2,usecols=(3,),unpack=True)
Flux_low  = np.loadtxt(data,skiprows=2,usecols=(8,),unpack=True)
Flux_up   = np.loadtxt(data,skiprows=2,usecols=(9,),unpack=True)

#/*****************************************************************************
# CURVE PLOTTING
#/*****************************************************************************
print 'Please wait...I am plotting for you'
plt.plot(E, (pro_spec*pow(E,2)), ls='-' , color='black',lw=2)
pam, = plt.plot(Emean_pam,Flux_pam,'ro')
plt.errorbar(Emean_pam,Flux_pam,xerr=[Emin_pam, Emax_pam],
 	         yerr=[Flux_low,Flux_up],fmt=None, ecolor='b')

plt.axis([5.e-2,5.e3, 1.e1,5.e3],interpolation='none')
plt.xlabel('E [GeV]',fontsize=20)
plt.ylabel(r'E$^2$ J(E)[GeV m$^{-2}$ s$^{-1}$ sr$^{-1}$]',fontsize=20)
plt.xscale('log')
plt.yscale('log')
plt.legend([pam],[r'\textbf{\textit{PAMELA 2009}}'],loc='upper right')
plt.savefig(OUTFOLD + '/' + 'pyprotons.eps',format='eps',bbox_inches='tight', dpi=100)
plt.show()

exit ()
