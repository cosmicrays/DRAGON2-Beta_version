#!/bin/bash/python

import numpy as np
import os
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

Run_Names = ['run_2D_breaks']
Number_of_Runs = len(Run_Names)
File_Names = []
for i in range(Number_of_Runs):
  File_Names.append("../output" + '/' + Run_Names[i] + '_spectrum.fits.gz')

norm = [1.0]                   # normalization factor
phi  = [0.2]                   # modulation potential
colors = ['black']

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

    B10     = np.zeros(dimE)
    B11     = np.zeros(dimE)
    C12     = np.zeros(dimE)
    N14     = np.zeros(dimE)
    O16     = np.zeros(dimE)


    B10_found = 0
    B11_found = 0
    C12_found = 0
    N14_found = 0
    O16_found = 0

    for i in xrange(1,n_ext):

        if hdulist[i].header['A'] == 10 and hdulist[i].header['Z_'] == 5 :
            B10_found = 1
            print "reading B10"
            current_HDU = hdulist[i].data
            for ie in xrange(0,dimE):
                B10[ie] = B10[ie] + current_HDU[ie]
        elif hdulist[i].header['A'] == 11 and hdulist[i].header['Z_'] == 5 :
            B11_found = 1
            print "reading B11"
            current_HDU = hdulist[i].data
            for ie in xrange(0,dimE):
                B11[ie] = B11[ie] + current_HDU[ie]
        elif hdulist[i].header['A'] == 12 and hdulist[i].header['Z_'] == 6 :
            C12_found = 1
            print "reading C12"
            current_HDU = hdulist[i].data
            for ie in xrange(0,dimE):
                C12[ie] = C12[ie] + current_HDU[ie]
        elif hdulist[i].header['A'] == 14 and hdulist[i].header['Z_'] == 7 :
            N14_found = 1
            print "reading N14"
            current_HDU = hdulist[i].data
            for ie in xrange(0,dimE):
                N14[ie] = N14[ie] + current_HDU[ie]
        elif hdulist[i].header['A'] == 16 and hdulist[i].header['Z_'] == 8 :
            O16_found = 1
            print "reading O16"
            current_HDU = hdulist[i].data
            for ie in xrange(0,dimE):
                O16[ie] = O16[ie] + current_HDU[ie]

    print 'Plotting...'

    B, = plt.loglog(E, E**2.8*(B10+B11), color="red", linestyle="-", label="Boron", lw=1.6)
    C, = plt.loglog(E, E**2.8*C12, color="blue", linestyle="-", label="Carbon", lw=1.6)
    N, = plt.loglog(E, E**2.8*N14, color="violet", linestyle="-", label="Nitrogen", lw=1.6)
    O, = plt.loglog(E, E**2.8*O16, color="navy", linestyle="-", label="Oxygen", lw=1.6)

plt.xlim(0.1,10000.)
plt.xlabel(r'E$_{k}$ [GeV/nuc]',fontsize=20)
plt.ylabel(r'Nuclei',fontsize=20)
plt.legend(loc='lower left', fontsize=16)

plt.savefig('nuclei.pdf',format='pdf',bbox_inches='tight', dpi=300)
plt.show()

exit ()

















