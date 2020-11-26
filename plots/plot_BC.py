#!/bin/bash/python

###############################################################################
# Licence: DRAGON TEAM                                                        #
# by G. Di Bernardo & D. Gaggero                                              # 
# Rev: September 2014                                                         #
###############################################################################

import os
import numpy as np
import math
from matplotlib import pyplot as plt
from matplotlib import rc, rcParams 
from scipy import interpolate
from scipy.interpolate import griddata  # useful for not-regular grid points
import pylab
import pyfits

# activate latex text rendering 
rc('text', usetex=True)
# Change all fonts to 'Computer Modern'
rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Sans Serif']})
rc('xtick', labelsize=18)
rc('ytick', labelsize=18)
rcParams['legend.numpoints'] = 1

print '##########################################################'
print 'Welcome to the B/C plotting routine                    '
print '##########################################################'

#/*****************************************************************************
# CONSTANTS: physical constants and conversion factors (mostly in cgs units)
#/*****************************************************************************

cmperpc = 3.08567758*1.e18
secperyr = 31556926.0
sec_yr = pow(secperyr,-1)       # to pass from seconds to year
pi = math.pi                    # pi-greek
pc_mt = cmperpc*1.e-2           # parsec in meters
m_to_pc= pow(pc_mt,-1)          # to pass from meters in parsec
m_p = .939                      # proton mass in GeV/c^2
m_e = 0.000510998918            # electron mass in GeV/c^2
b = 1.4e-16                     # energy loss rate in GeV^-1*s^-1
ergperev = 1.6021773*1.e-12
GeV_erg = ergperev*1.e9         # from GeV to erg
c2kB = 6.50966198               # divided by 10^36

#/*****************************************************************************
# FILENAMES
#/*****************************************************************************

Input_Folder    = os.path.join(os.environ['HOME'],'python','runs')
Output_Folder   = os.path.join(os.environ['HOME'],'python','plots')
Data_Folder     = os.path.join(os.environ['HOME'],'python','data')

Run_Names = ['run_2D', 'run_2D', 'run_2D']
Number_of_Runs = len(Run_Names)
File_Names = []
for i in range(Number_of_Runs):
  File_Names.append(Input_Folder + '/' + Run_Names[i] + '_spectrum.fits.gz')

norm = [1.0, 1.0, 1.0]                   # normalization factor
phi  = [0.2, 0.4, 0.6]                   # modulation potential
colors = ['black', 'red', 'blue']


def hadron_modulation(E, Elog, spectrum, potential, A_, Z_):
	print "modulating nucleus with A=", A_, " Z=",Z_
	print "potential=",potential
	dimE=len(E)
	spectrum_mod =  np.zeros(int(dimE))
	A = ((E + m_p)**2 - m_p**2) / (((E + m_p) + (A_/Z_)*potential)**2 - m_p**2)
	F = interpolate.InterpolatedUnivariateSpline(Elog, spectrum)
	Elog_shifted = np.zeros(dimE)
	shift = 0.
	shift = (A_/Z_)*potential
	for i in xrange(dimE):
	  	Elog_shifted[i] = np.log(E[i] + shift)
	spectrum_mod = A * F(Elog_shifted)
	return spectrum_mod	


#/*****************************************************************************
# COMPUTATION
#/*****************************************************************************

for k in range(Number_of_Runs):

	print " "
	print "Reading file ", File_Names[k]
	print " "

	hdulist = pyfits.open(File_Names[k]) 
	n_ext  = len(hdulist)

	table_hdu    = hdulist[0]              
	table_header = table_hdu.header     # header attribute of TABLE 
	emin   = table_header['EKMIN']
	ek_fac = table_header['EKIN_FAC']
	dimE   = table_header['DIME']
	energy = np.zeros(dimE,float)

	E = np.zeros(dimE)
	for i in xrange(0,dimE):
		E[i] = emin*pow(ek_fac,i)

	Elog = np.zeros(dimE)
	for i in xrange(0,dimE):
		Elog[i] = np.log(E[i])

	table_B10     = np.zeros(dimE)
	table_B11     = np.zeros(dimE)
	table_C12     = np.zeros(dimE)
	table_C13     = np.zeros(dimE)
	table_C14     = np.zeros(dimE)
	B10     = np.zeros(dimE)
	B11     = np.zeros(dimE)
	C12     = np.zeros(dimE)
	C13     = np.zeros(dimE)
	C14     = np.zeros(dimE)
	B10_mod     = np.zeros(dimE)
	B11_mod     = np.zeros(dimE)
	C12_mod     = np.zeros(dimE)
	C13_mod     = np.zeros(dimE)
	C14_mod     = np.zeros(dimE)

        B10_found = 0
        B11_found = 0
        C12_found = 0
	C13_found = 0
	C14_found = 0

	for i in xrange(1,n_ext):

		if hdulist[i].header['A'] == 10 and hdulist[i].header['Z_'] == 5 :
			B10_found = 1
			print "reading B10"
			current_HDU = hdulist[i].data
			for ie in xrange(0,dimE): 
				table_B10[ie] = table_B10[ie] + current_HDU[ie]
		elif hdulist[i].header['A'] == 11 and hdulist[i].header['Z_'] == 5 :
			B11_found = 1
			print "reading B11"
			current_HDU = hdulist[i].data
			for ie in xrange(0,dimE): 
				table_B11[ie] = table_B11[ie] + current_HDU[ie]
		elif hdulist[i].header['A'] == 12 and hdulist[i].header['Z_'] == 6 :
			C12_found = 1
			print "reading C12"
			current_HDU = hdulist[i].data
			for ie in xrange(0,dimE): 
				table_C12[ie] = table_C12[ie] + current_HDU[ie]
		elif hdulist[i].header['A'] == 13 and hdulist[i].header['Z_'] == 6 :
			C13_found = 1
			print "reading C13"
			current_HDU = hdulist[i].data
			for ie in xrange(0,dimE): 
				table_C13[ie] = table_C13[ie] + current_HDU[ie]
		elif hdulist[i].header['A'] == 14 and hdulist[i].header['Z_'] == 6 :
			C14_found = 1
			print "reading C14"
			current_HDU = hdulist[i].data
			for ie in xrange(0,dimE): 
				table_C14[ie] = table_C14[ie] + current_HDU[ie]	

	if (B10_found == 1):
		print "I am modulating boron 10..."
		B10 = table_B10
		B10_mod = hadron_modulation(E, Elog, B10, phi[k], 5., 10.)
	else:
		print "Warning: B10 not found in ", Run_Names[k]

	if (B11_found == 1):
		print "I am modulating boron 11..."
		B11 = table_B11
		B11_mod = hadron_modulation(E, Elog, B11, phi[k], 5., 11.)
	else:
		print "Warning: B11 not found in ", Run_Names[k]

	if (C12_found == 1):			
		print "I am modulating carbon 12..."
		C12 = table_C12
		C12_mod = hadron_modulation(E, Elog, C12, phi[k], 6., 12.)
	else:
		print "Warning: C12 not found in ", Run_Names[k]

	if (C13_found == 1):			
		print "I am modulating carbon 13..."
		C13 = table_C13
		C13_mod = hadron_modulation(E, Elog, C13, phi[k], 6., 13.)
	else:
		print "Warning: C13 not found in ", Run_Names[k]
	
	if (C14_found == 1):	
		print "I am modulating carbon 14..."
		C14 = table_C14
		C14_mod = hadron_modulation(E, Elog, C14, phi[k], 6., 14.)
	else:	
		print "Warning: C14 not found in ", Run_Names[k]

	boron  = (B10 + B11)
	carbon = (C12 + C13 + C14)
	BoverC = boron/carbon 
	B_mod   = (B10_mod + B11_mod)           # modulated boron spectrum 
	C_mod   = (C12_mod + C13_mod + C14_mod) # modulated carbon spectrum 
	BoverC_mod = (B_mod/C_mod)               # modulated boron to carbon ratio	

	#/*****************************************************************************
	# PLOTTING
	#/*****************************************************************************

	print 'Plotting...'

	BonC_mod, = plt.plot(E, BoverC_mod, ls='-' ,  color=colors[k],lw=2)
	BonC,     = plt.plot(E, BoverC, ls='--' , color=colors[k],lw=2)


#/*****************************************************************************
# HANDLING EXPERIMENTAL DATA
#/*****************************************************************************

data = Data_Folder + '/' + 'bc_ams2.dat' # preliminary AMS data 
				      # Emin,Emax = np.loadtxt(data,skiprows=2,usecols=(1,2),unpack=True)
Emean_ams     = np.loadtxt(data,skiprows=1,usecols=(0,),unpack=True)
Ratio_ams     = np.loadtxt(data,skiprows=1,usecols=(2,),unpack=True)
Ratio_low_ams = np.loadtxt(data,skiprows=1,usecols=(1,),unpack=True)
Ratio_up_ams  = np.loadtxt(data,skiprows=1,usecols=(3,),unpack=True)

data = Data_Folder + '/' + 'PAMELA_BC.dat' # preliminary AMS data 
				      # Emin,Emax = np.loadtxt(data,skiprows=2,usecols=(1,2),unpack=True)
Emean_pam     = np.loadtxt(data,skiprows=1,usecols=(0,),unpack=True)
Ratio_pam     = np.loadtxt(data,skiprows=1,usecols=(3,),unpack=True)
Err_low_pam = np.loadtxt(data,skiprows=1,usecols=(8,),unpack=True)
Err_up_pam  = np.loadtxt(data,skiprows=1,usecols=(9,),unpack=True)
Ratio_low_pam = Ratio_pam - Err_low_pam
Ratio_up_pam  = Ratio_pam + Err_up_pam

#/*****************************************************************************
#  FINALIZING PLOT
#/*****************************************************************************

print " "
print "Finalizing plot..."
print " "

ams02, = plt.plot(Emean_ams,Ratio_ams,'o',color='grey')
plt.errorbar(Emean_ams,Ratio_ams,
   	         yerr=[Ratio_ams-Ratio_low_ams,Ratio_up_ams-Ratio_ams],fmt=None, ecolor='grey')

pam, = plt.plot(Emean_pam,Ratio_pam,'v',color='orange',ms=8)
plt.errorbar(Emean_pam,Ratio_pam,
   	         yerr=[Ratio_pam-Ratio_low_pam,Ratio_up_pam-Ratio_pam],fmt=None, ecolor='orange')
 
plt.axis([3.e-2,3.e3, 0.,0.4],interpolation='none')
plt.xlabel(r'E$_{k}$ [GeV/nuc]',fontsize=20)
plt.ylabel(r'B/C',fontsize=20)
plt.xscale('log')

#legend1 = plt.legend([BonC,BonC2,BonC3],[r'\textbf{KRA model}',r'\textbf{PD model}',r'\textbf{KOL model}'])
legend2 = plt.legend([ams02,pam],
			[r'\textbf{AMS-02}',r'\textbf{PAMELA}'],
			#bbox_to_anchor = (1.5,1.),
			loc='lower left')
#legend2 = plt.legend([BonC],[r'$D_{0} = 3.0$, $\chi^{2}_{AMS} = %f$'% chi2_ams],
#			loc='lower left')
#plt.gca().add_artist(legend1)
plt.savefig(Output_Folder + '/'  + 'BC.eps',format='eps',bbox_inches='tight', dpi=100)
plt.show()

exit ()

















