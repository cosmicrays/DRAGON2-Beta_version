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
#from scipy.interpolate import interp1d  
import pylab
import pyfits

# activate latex text rendering 
rc('text', usetex=True)
# Change all fonts to 'Computer Modern'
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('xtick', labelsize=18)
rc('ytick', labelsize=18)
rcParams['legend.numpoints'] = 1

print '##########################################################'
print 'Welcome to the proton plotting routine                    '
print '##########################################################'

#/*****************************************************************************
# CONSTANTS: physical constants and conversion factors (mostly in cgs units)
#/*****************************************************************************

pi = math.pi                    # pi-greek
m_p = .939                      # proton mass in GeV/c^2
m_e = 0.000510998918            # electron mass in GeV/c^2
b = 1.4e-16                     # energy loss rate in GeV^-1*s^-1
c2kB = 6.50966198               # divided by 10^36

#/*****************************************************************************
# FILENAMES
#/*****************************************************************************

Input_Folder    = os.path.join(os.environ['HOME'],'python','runs')
Output_Folder   = os.path.join(os.environ['HOME'],'python','plots')
Data_Folder     = os.path.join(os.environ['HOME'],'python','data')

Run_Names = ['run_2D','run_2D','run_2D']
Number_of_Runs = len(Run_Names)
File_Names = []
for i in range(Number_of_Runs):
  File_Names.append(Input_Folder + '/' + Run_Names[i] + '_spectrum.fits.gz')

norm = [1.0,1.0,1.0]                 # normalization factor
phi  = [0.2,0.4,0.6]                   # modulation potential
colors = ['red','green','blue']

for k in range(Number_of_Runs):

	print " "
	print "Reading file ", File_Names[k]
	print " "

	#/*****************************************************************************
	# Reading and modulating data
	#/*****************************************************************************

	hdulist = pyfits.open(File_Names[k]) 
	n_ext  = len(hdulist)

	# Header keywords

	table_hdu    = hdulist[0]              
	table_header = table_hdu.header     
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
	 
	# Energy scale in GeV 

	E = np.zeros(dimE)
	for i in xrange(0,dimE):
		E[i] = emin*pow(ek_fac,i)

	Elog = np.zeros(dimE)
	for i in xrange(0,dimE):
		Elog[i] = np.log(E[i])

	# Reading protons
	
	print "looking for protons..."
	table_protons = np.zeros(dimE)
	for i in xrange(1,n_ext):
	        print "A = ",hdulist[i].header['A']," Z = ", hdulist[i].header['Z_']," is secondary = ", hdulist[i].header['SEC']
		if hdulist[i].header['A'] == 1 and hdulist[i].header['Z_'] == 1:
			if hdulist[i].header['SEC'] == 0 :
				print "These are primary protons!"
			if hdulist[i].header['SEC'] == 1 :
				print "These are secondary protons!"
			current_HDU = hdulist[i].data 
			for ie in xrange(0,dimE):
				table_protons[ie] = table_protons[ie] + current_HDU[ie]


	print "Modulating protons..."

	protons = np.zeros(dimE)
	protons_mod = np.zeros(dimE) 

	protons     = norm[k]*table_protons

	A = ((E + m_p)**2 - m_p**2) / (((E + m_p) + phi[k])**2 - m_p**2)
	
	F = interpolate.InterpolatedUnivariateSpline(Elog, protons)
	Elog_shifted = np.zeros(dimE)
	for i in xrange(dimE):
	  	Elog_shifted[i] = np.log(E[i] + phi[k])
	protons_mod = A * F(Elog_shifted)

	print 'Modulated protons: ',protons_mod 

	#/*****************************************************************************
	# PLOTTING
	#/*****************************************************************************

	print 'Plotting...'

	spec, = plt.plot(E, (protons_mod*pow(E,2)), ls='-' , color=colors[k],lw=2)
	lis,  = plt.plot(E, (protons*pow(E,2)), ls='--' , color=colors[k],lw=2)

#/*****************************************************************************
# READING EXPERIMENTAL DATA
#/*****************************************************************************

data = Data_Folder + '/' + 'PAMELA_protons.dat'
Emin_pam,Emax_pam = np.loadtxt(data,skiprows=2,usecols=(1,2),unpack=True)
Emean_pam = np.loadtxt(data,skiprows=2,usecols=(0,),unpack=True)
Flux_pam  = np.loadtxt(data,skiprows=2,usecols=(3,),unpack=True)
#Flux_low  = np.loadtxt(data,skiprows=2,usecols=(8,),unpack=True)
#Flux_up   = np.loadtxt(data,skiprows=2,usecols=(9,),unpack=True)
Err_low_pam = np.loadtxt(data,skiprows=1,usecols=(8,),unpack=True)
Err_up_pam  = np.loadtxt(data,skiprows=1,usecols=(9,),unpack=True)
Flux_low_pam = Flux_pam - Err_low_pam
Flux_up_pam  = Flux_pam + Err_up_pam
 
E2xFlux_pam = ((Emean_pam)**2) * Flux_pam
Emin_pam_err = (Emean_pam - Emin_pam)
Emax_pam_err = (Emax_pam - Emean_pam)

#data = Data_Folder + '/' + 'protons_AMS02.dat'
##Emin_pam,Emax_pam = np.loadtxt(data,skiprows=1,usecols=(1,2),unpack=True)
#Rmean_ams = np.loadtxt(data,skiprows=2,usecols=(0,),unpack=True)
#Flux_ams  = np.loadtxt(data,skiprows=2,usecols=(2,),unpack=True)
#Flux_low  = np.loadtxt(data,skiprows=2,usecols=(1,),unpack=True)
#Flux_up   = np.loadtxt(data,skiprows=2,usecols=(3,),unpack=True)

#/*****************************************************************************
#  FINALIZING PLOT
#/*****************************************************************************

print " "
print "Finalizing plot..."
print " "

#pam, = plt.plot(Emean_pam,Flux_pam,'ro')
#plt.errorbar(Emean_pam,Flux_pam,xerr=[Emin_pam, Emax_pam],
#     yerr=[Flux_low,Flux_up],fmt=None, ecolor='b')
#ams, = plt.plot(Rmean_ams,Flux_ams,'ro')
#plt.errorbar(Rmean_ams,Flux_ams, yerr=[Flux_ams-Flux_low,Flux_up-Flux_ams],
#	     fmt=None, ecolor='b')

pam, = plt.plot(Emean_pam,E2xFlux_pam,'ro')
plt.errorbar(Emean_pam,E2xFlux_pam,xerr=[Emin_pam_err, Emax_pam_err],
     yerr=[Err_low_pam*pow(Emean_pam,2),Err_up_pam*pow(Emean_pam,2)],fmt=None, ecolor='b')
#ams, = plt.plot(Rmean_ams,Flux_ams,'ro')
#plt.errorbar(Rmean_ams,Flux_ams, yerr=[Flux_ams-Flux_low,Flux_up-Flux_ams],
#        fmt=None, ecolor='b')
 

plt.axis([1.e-1,5.e3, 1.e1,5.e3],interpolation='none')
#plt.axis([5.e-2,5.e3, 1.e-5,3.e4],interpolation='none')
plt.xlabel('E [GeV]',fontsize=20)
plt.ylabel(r'E$^2$ J(E)[GeV m$^{-2}$ s$^{-1}$ sr$^{-1}$]',fontsize=20)
#plt.ylabel(r'J(E)[GeV$^{-1}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]',fontsize=20)
plt.xscale('log')
plt.yscale('log')

plt.legend([pam],[r'\textbf{\textit{PAMELA 2009}}'],loc='upper right')
#plt.legend([spec],['$\Phi = 0.30$ GV'],loc='lower right')
#plt.legend([ams],[r'\textbf{\textit{AMS 02}}'],loc='upper right')

plt.savefig(Output_Folder + '/' + 'protons.eps',format='eps',bbox_inches='tight', dpi=200)
plt.show()

exit ()
