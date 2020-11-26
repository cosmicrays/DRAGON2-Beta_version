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
print 'Welcome to the lepton plotting routine                    '
print '##########################################################'

#/*****************************************************************************
# CONSTANTS: physical constants and conversion factors (mostly in cgs units)
#/*****************************************************************************

pi = math.pi                    # pi-greek
pc_mt = 3.0857e+16           # parsec in meters
m_to_pc= pow(pc_mt,-1)          # to pass from meters in parsec
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

Run_Names = ['run_2D', 'run_2D', 'run_2D']
Number_of_Runs = len(Run_Names)
File_Names = []
for i in range(Number_of_Runs):
  File_Names.append(Input_Folder + '/' + Run_Names[i] + '_spectrum.fits.gz')

norm = [1.0, 1.0, 1.0]                   # normalization factor
norm_extra = [1.0, 1.0, 1.0]             # normalization factor for the extra component   
phi  = [0.5, 0.6, 0.7]                # modulation potential
colors_positrons = ['black', 'red', 'blue']
colors_electrons = ['black', 'red', 'blue']
colors_extra_component = ['black', 'red', 'blue']

def lepton_modulation(E, Elog, spectrum, potential):
	dimE=len(E)
	spectrum_mod =  np.zeros(int(dimE))
	A = ((E + m_e)**2 - m_e**2) / (((E + m_e) + potential)**2 - m_e**2)
	F = interpolate.InterpolatedUnivariateSpline(Elog, spectrum)
	Elog_shifted = np.zeros(dimE)
	for i in xrange(dimE):
	  	Elog_shifted[i] = np.log(E[i] + potential)
	spectrum_mod = A * F(Elog_shifted)
	return spectrum_mod	

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

	table_pos     = np.zeros(dimE)
	table_pri_ele = np.zeros(dimE)
	table_sec_ele = np.zeros(dimE)
	table_extra   = np.zeros(dimE)

	sec_pos_found = 0
	sec_ele_found = 0
	pri_ele_found = 0
	extra_comp_found = 0

	for i in xrange(1,n_ext):

		if hdulist[i].header['A'] == 0 and hdulist[i].header['Z_'] == 1  and hdulist[i].header['EXTRA'] == 0 and hdulist[i].header['SEC'] == 1 and hdulist[i].header['DM'] == 0:
			sec_pos_found = 1
			print 'reading secondary positrons'
			current_HDU = hdulist[i].data 
			for ie in xrange(0,dimE):
				table_pos[ie] = table_pos[ie] + current_HDU[ie]
#
		elif hdulist[i].header['A'] == 0 and hdulist[i].header['Z_'] == -1 and hdulist[i].header['EXTRA'] == 0 and hdulist[i].header['SEC'] == 0 and hdulist[i].header['DM'] == 0:
			pri_ele_found = 1
			print 'reading primary electrons'
			current_HDU = hdulist[i].data 
			for ie in xrange(0,dimE):
				table_pri_ele[ie] = table_pri_ele[ie] + current_HDU[ie]	
#
		elif hdulist[i].header['A'] == 0 and hdulist[i].header['Z_'] == -1 and hdulist[i].header['EXTRA'] == 0 and hdulist[i].header['SEC'] == 1 and hdulist[i].header['DM'] == 0:
			sec_ele_found = 1
			print 'reading secondary electrons'
			current_HDU = hdulist[i].data 
			for ie in xrange(0,dimE):
				table_sec_ele[ie] = table_sec_ele[ie] + current_HDU[ie]	
#
		elif hdulist[i].header['A'] == 0 and hdulist[i].header['Z_'] == -1 and hdulist[i].header['EXTRA'] == 1 and hdulist[i].header['SEC'] == 0 and hdulist[i].header['DM'] == 0:
			extra_comp_found = 1
			print 'reading extra component'
			current_HDU = hdulist[i].data 
			for ie in xrange(0,dimE):
				table_extra[ie] = table_extra[ie] + current_HDU[ie]	

	if (sec_pos_found == 0):
		print "Warning: secondary positrons not found in ", Run_Names[k]
	if (sec_ele_found == 0):
		print "Warning: secondary electrons not found in ", Run_Names[k]
	if (pri_ele_found == 0):
		print "Warning: primary electrons not found in ", Run_Names[k]
	if (extra_comp_found == 0):
		print "Warning: extra component not found in ", Run_Names[k]

	total_positrons = np.zeros(dimE)
	total_positrons_mod = np.zeros(dimE)
	for ie in xrange(0,dimE):
		total_positrons[ie] = total_positrons[ie] + table_pos[ie] + norm_extra[k]*table_extra[ie] 
	total_electrons = np.zeros(dimE)
	total_electrons_mod = np.zeros(dimE)
	for ie in xrange(0,dimE):
		total_electrons[ie] = total_electrons[ie] + norm[k]*table_pri_ele[ie] + table_sec_ele[ie] + norm_extra[k]*table_extra[ie] 
	extra_component = np.zeros(dimE)	
	for ie in xrange(0,dimE):
		extra_component[ie] = norm_extra[k]*table_extra[ie]

	print "Modulating total positrons..."
	total_positrons_mod = lepton_modulation(E, Elog, total_positrons, phi[k])
	print "Modulating total electrons..."
	total_electrons_mod = lepton_modulation(E, Elog, total_electrons, phi[k])

	#/*****************************************************************************
	# PLOTTING
	#/*****************************************************************************

	print 'Plotting...'

	spec_pos, = plt.plot(E, (total_positrons_mod*pow(E,3)), ls='-' , color=colors_positrons[k],lw=2)
	lis_pos,  = plt.plot(E, (total_positrons*pow(E,3)), ls='--' , color=colors_positrons[k],lw=2)
	spec_ele, = plt.plot(E, (total_electrons_mod*pow(E,3)), ls='-' , color=colors_electrons[k],lw=2)
	lis_ele,  = plt.plot(E, (total_electrons*pow(E,3)), ls='--' , color=colors_electrons[k],lw=2)
	lis_extra,  = plt.plot(E, (extra_component*pow(E,3)), ls='-.' , color=colors_extra_component[k],lw=2)


#/*****************************************************************************
# HANDLING EXPERIMENTAL DATA
#/*****************************************************************************

data = Data_Folder + '/' + 'AMS02_PFraction.dat'
Emin_ams,Emax_ams = np.loadtxt(data,skiprows=4,usecols=(0,1),unpack=True)
Ratio_ams = np.loadtxt(data,skiprows=4,usecols=(3,),unpack=True)
Ratio_stat, Ratio_sys = np.loadtxt(data,skiprows=4,usecols=(4,10),unpack=True)
Emean_ams = (Emax_ams + Emin_ams) / 2.
Emin_ams_err = (Emean_ams - Emin_ams)
Emax_ams_err = (Emax_ams - Emean_ams)
Tot_ams_err  = np.sqrt(Ratio_sys**2 + Ratio_stat**2)

data = Data_Folder + '/' + 'AllElectronFlux_AMS02_fromPlots.dat'
#Emin_pam,Emax_pam = np.loadtxt(data,skiprows=1,usecols=(1,2),unpack=True)
Emean_e_ams02     = np.loadtxt(data,skiprows=2,usecols=(0,),unpack=True)
Flux_e_ams02      = np.loadtxt(data,skiprows=2,usecols=(2,),unpack=True)
Flux_e_low_ams02  = np.loadtxt(data,skiprows=2,usecols=(1,),unpack=True)
Flux_e_up_ams02   = np.loadtxt(data,skiprows=2,usecols=(3,),unpack=True)

data = Data_Folder + '/' + 'PositronFlux_AMS02_fromPlots.dat'
#Emin_pam,Emax_pam = np.loadtxt(data,skiprows=1,usecols=(1,2),unpack=True)
Emean_p_ams02     = np.loadtxt(data,skiprows=2,usecols=(0,),unpack=True)
Flux_p_ams02      = np.loadtxt(data,skiprows=2,usecols=(2,),unpack=True)
Flux_p_low_ams02  = np.loadtxt(data,skiprows=2,usecols=(1,),unpack=True)
Flux_p_up_ams02   = np.loadtxt(data,skiprows=2,usecols=(3,),unpack=True)

data = Data_Folder + '/' + 'LEPTONS_DATAFILE.dat'
ID_exp = np.loadtxt(data,skiprows=30,usecols=(0,),unpack=True)
Emin,Emax = np.loadtxt(data,skiprows=30,usecols=(2,3),unpack=True)
Emean,Flux = np.loadtxt(data,skiprows=30,usecols=(1,4),unpack=True)
Flux_stat,Flux_sys = np.loadtxt(data,skiprows=30,usecols=(5,7),unpack=True)
ID_pam13 = np.where(ID_exp==11)
Emin_pam13 = Emin[ID_pam13]
Emax_pam13 = Emax[ID_pam13]
Emean_pam13 = Emean[ID_pam13]
Flux_pam13 = Flux[ID_pam13]
Flux_pstat13 = Flux_stat[ID_pam13]
Flux_psys13 = Flux_sys[ID_pam13]
Emin_pam13_err = (Emean_pam13 - Emin_pam13)
Emax_pam13_err = (Emax_pam13 - Emean_pam13)
E3xFlux_pam13 = ((Emean_pam13)**3) * Flux_pam13
delta_F = np.sqrt(Flux_pstat13**2 + Flux_psys13**2)
# in the propagation error, we assume delta_E = 0 
tot_err = delta_F*(Emean_pam13**3)

ID_pam13 = np.where(ID_exp==12) # this time we select the positron ratio data
Em_rpam13 = Emean[ID_pam13]
Ratio_pam13 = Flux[ID_pam13]
Ratio_pstat13 = Flux_stat[ID_pam13]

#/*****************************************************************************
#  FINALIZING PLOT
#/*****************************************************************************

print " "
print "Finalizing plot..."
print " "

ams02_e, = plt.plot(Emean_e_ams02,Flux_e_ams02,'ro', color='grey')
plt.errorbar(Emean_e_ams02,Flux_e_ams02,
			yerr=[Flux_e_ams02 - Flux_e_low_ams02,Flux_e_up_ams02 - Flux_e_ams02],
 	        ecolor='k',fmt=None)

ams02_p, = plt.plot(Emean_p_ams02,Flux_p_ams02,'o', color='grey')
plt.errorbar(Emean_p_ams02,Flux_p_ams02,
			yerr=[Flux_p_ams02 - Flux_p_low_ams02,Flux_p_up_ams02 - Flux_p_ams02],
 	        ecolor='k',fmt=None)

pam13, = plt.plot(Emean_pam13,E3xFlux_pam13,'ro',color='red', label=r'PAMELA 2013 [e$^{+}$]')
plt.errorbar(Emean_pam13,E3xFlux_pam13,xerr=[Emin_pam13_err, Emax_pam13_err],
 	         yerr=tot_err,fmt=None, ecolor='r')


plt.axis([5.e-1,1.e4, 1.e0,3.e3],interpolation='none')
plt.xlabel('E [GeV]',fontsize=20)
plt.ylabel('E$^3$ J(E)[GeV$^2$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]',fontsize=20)
plt.xscale('log')
plt.yscale('log')

plt.legend([pam13,ams02_e,ams02_p],[
 	       r'\textbf{PAMELA 2013 [e$^{+}$]}',r'\textbf{AMS-02 2013 [e$^{-}$]}',
 	       r'\textbf{AMS-02 2013 [e$^{+}$]}'],loc='upper right')

plt.savefig(Output_Folder + '/' + 'leptons.eps',format='eps',bbox_inches='tight', dpi=200)
plt.show()
plt.hold('off')

exit()

