#!/bin/bash/python

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

Input_Folder    = "/Users/daniele/Dropbox/Fisica/2018/SFRproject/SF_GC/DRAGON-SF/output/"
Output_Folder   = "/Users/daniele/Dropbox/Fisica/2018/SFRproject/SF_GC/DRAGON-SF/plots/"
Data_Folder     = "/Users/daniele/Dropbox/Fisica/2018/SFRproject/SF_GC/DRAGON-SF/plots/"

Run_Names = ['run_3D_spike']
Number_of_Runs = len(Run_Names)
File_Names = []
for i in range(Number_of_Runs):
  File_Names.append(Input_Folder + '/' + Run_Names[i] + '.fits.gz')

norm = [1.0]                 # normalization factor
colors = ['black', 'red', 'blue']

# which particle do you want to plot?
A_ = 0
Z_ = -1

# which plane do you want to plot?
z_planes = [0.]

# which vertical profiles do you want to plot?
R_positions = [0.5]

# which horizontal profiles do you want to plot?
z_positions = [0.]

# which energy do you want to plot?
E_plot = [1.]

list_of_contours = []
list_of_x_vectors = []
list_of_y_vectors = []
list_of_z_vectors = []
list_of_z_profiles = []
list_of_R_profiles = []

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
	xmin   = table_header['xmin']
	xmax   = table_header['xmax']
	ymin   = table_header['ymin']
	ymax   = table_header['ymax']
	zmin   = table_header['zmin']
	zmax   = table_header['zmax']
	ixsun  = table_header['ixsun']
	iysun  = table_header['iysun']
	izsun  = table_header['izsun']

	print "ixsun = ", ixsun
	print "iysun = ", iysun
	print "izsun = ", izsun

	dimx   = table_header['dimx']
	dimy   = table_header['dimy']
	dimz   = table_header['dimz']
	xobs   = table_header['xobs']
	yobs   = table_header['yobs']
	zobs   = table_header['zobs']
	emin   = table_header['ekmin']
	ek_fac = table_header['ekin_fac']
	dimE   = table_header['dimE'] 
	 
	# Energy scale in GeV 

	print "Energy scale"
	E = np.zeros(dimE)
	for i in xrange(0,dimE):
		E[i] = emin*pow(ek_fac,i)
		print E[i]

	Elog = np.zeros(dimE)
	for i in xrange(0,dimE):
		Elog[i] = np.log(E[i])

	ie = 0
	while E[ie] < E_plot[k]:
		ie += 1

	print "Current energy = ", E[ie]

	dx = (xmax - xmin) / (dimx - 1)
	dy = (ymax - ymin) / (dimy - 1)
	dz = (zmax - zmin) / (dimz - 1)
	x = xmin + np.arange(dimx)*dx
	y = ymin + np.arange(dimy)*dy
	z = zmin + np.arange(dimz)*dz

	list_of_x_vectors.append(x)
	list_of_y_vectors.append(y)
	list_of_z_vectors.append(z)

	iz = 0
	while z[iz] <= z_planes[k]:
		iz +=1

	ix_profile = 0
	while x[ix_profile] <= R_positions[k]:
		ix_profile +=1

	iz_profile = 0
	while z[iz_profile] <= z_positions[k]:
		iz_profile +=1

	# Reading the particle density from FITS file
	
	print "Looking for the particle with A=", A_, " and Z=", Z_
	particle_density = np.zeros((dimz,dimy,dimx,dimE))

	for i in xrange(1,n_ext):

		print "A = ",hdulist[i].header['A']," Z = ", hdulist[i].header['Z_']," is secondary = ", hdulist[i].header['SEC']
		if hdulist[i].header['A'] == A_ and hdulist[i].header['Z_'] == Z_:

				print "Particle found!"
				current_HDU = hdulist[i].data
				print "Reading HDU with dimension ", current_HDU.shape
				particle_density = particle_density + current_HDU

	print "Datacube"
	print particle_density

	current_map_2D = (particle_density[iz,:,:,ie]) / np.max(particle_density[iz,:,:,ie])
	list_of_contours.append(current_map_2D)

	current_R_profile = (particle_density[iz_profile,iysun,:,ie]) #/ np.max(particle_density[iz_profile,iysun,:,ie])

	print "R_profile"
	print current_R_profile

	list_of_R_profiles.append(current_R_profile)

	current_z_profile = (particle_density[:,iysun,ix_profile,ie]) #/ np.max(particle_density[:,iysun,ix_profile,ie])
	list_of_z_profiles.append(current_z_profile)


#/*****************************************************************************
# PLOTTING
#/*****************************************************************************

"""
print " "
print "Plotting contours..."
print " "

for k in range(Number_of_Runs):

	plt.figure()
	levels = 50
	X, Y = np.meshgrid(list_of_x_vectors[k], list_of_y_vectors[k])
	plot1 =  plt.contour(X, Y, list_of_contours[k], levels, colors='black', linewidth = .5)
	plot2 = plt.contourf(X, Y, list_of_contours[k], alpha = 0.75, cmap = 'Spectral')

	plt.xlabel('X [kpc]',fontsize=16)
	plt.ylabel('Y [kpc]',fontsize=16)
	cbar = plt.colorbar(plot2)
	cbar.set_label(r'Proton density, E$_{p}$ = 1 GeV, z = 0', fontsize=16)
	plt.savefig(Output_Folder + '/' + 'contour_' + str(k) + '.eps',format='eps',dpi=300)
	plt.show()
"""


print " "
print "Plotting radial profiles..."
print " "

plt.figure()
max_x=-1000.
for k in range(Number_of_Runs):	
	profile, = plt.plot(list_of_x_vectors[k],list_of_R_profiles[k],ls='-' , color=colors[k],lw=2)
	current_max_x=np.max(list_of_x_vectors[k])
	if (current_max_x > max_x):
		max_x = current_max_x
#plt.axis([-max_x,max_x, 0.,1.5],interpolation='none')
plt.xlabel('x [kpc]',fontsize=20)
plt.ylabel('profile',fontsize=20)
plt.savefig(Output_Folder + '/' + 'radial_profile_' + str(k) + '.pdf',format='pdf',dpi=300)
plt.show()
		
print " "
print "Plotting vertical profiles..."
print " "

plt.figure()
min_z=1000.
max_z=-1000.
for k in range(Number_of_Runs):
	profile, = plt.plot(list_of_z_vectors[k],list_of_z_profiles[k],ls='-' , color=colors[k],lw=2)
	current_min_z = np.min(list_of_z_vectors[k])
	current_max_z = np.max(list_of_z_vectors[k])
	if (current_min_z < min_z):
		min_z = current_min_z
	if (current_max_z > max_z):
		max_z = current_max_z
#plt.axis([min_z,max_z, 0.,1.5],interpolation='none')
plt.xlabel('z [kpc]',fontsize=20)
plt.ylabel('profile',fontsize=20)
plt.savefig(Output_Folder + '/' + 'vertical_profile_' + str(k) + '.pdf',format='pdf',dpi=300)
plt.show()
	

exit ()
