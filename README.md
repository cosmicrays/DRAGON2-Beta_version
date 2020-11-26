# DRAGON2 - BETA VERSION
Diffusion Reacceleration and Advection of Galactic cosmic rays: an Open New code
Version 2 Beta


## References

We refer to the following technical papers:

[PAPER I] Carmelo Evoli, Daniele Gaggero, Andrea Vittino, Giuseppe Di Bernardo, Mattia Di Mauro, Arianna Ligorini, Piero Ullio, Dario Grasso, "Cosmic-ray propagation with DRAGON2: I. numerical solver and astrophysical ingredients", JCAP 02 (2017) 015 (https://arxiv.org/abs/1607.07886)

[PAPER II] Carmelo Evoli, Daniele Gaggero, Andrea Vittino, Mattia Di Mauro, Dario Grasso, Mario Nicola Mazziotta, "Cosmic-ray propagation with DRAGON2: II. Nuclear interactions with the interstellar gas", JCAP 07 (2018) 006  (https://arxiv.org/abs/1711.09616)


## Installation

## 1) Setting up the system

External LIBRARIES needed:  

- [GSL] (http://www.gnu.org/software/gsl/) 
- [CFITSIO] (http://heasarc.gsfc.nasa.gov/fitsio/) 

Please make sure you have the GNU autotools and the GNU gcc compiler properly installed.

HINT for MAC users:  We recommend to install the GCC compiler, the GNU autotools and the GSL/cfitsio libraries consistently with the same package manager. 
In particular, the installation process is fully tested within the HomeBrew environment ( https://brew.sh/ )
Within Homebrew, the following packages may be needed:
- autoconf
- automake
- libtool
- cfitsio
- gsl
- gcc

## 2) Initialization

Before installing the code launch the script to initialize installation tools:

`./start.sh

For MAC users:

./startMAC.sh`

## 3) configure

Configure the code, a typical command line is:

`./configure --with-cfitsio=$CFITSIO_DIR --with-numcpu=2`

HINT for MAC users: We recommend to use the GCC compiler (e.g. as provided by HomeBrew) and explicitly instruct configure to use it. Please check this example:
./configure --with-cfitsio=$CFITSIO_DIR  CXX=g++-9 CC=gcc-9 --with-numcpu=4
 
where `$CFITSIO_DIR` is the path of your cfitsio library and `NUMCPU` is the machine core number.

The default installation path is in the same folder as the source code is (the program automatically creates the `bin/` and `lib/` subfolders). It can be set via `--prefix=<NEW_INSTALLATION_PATH>`

## 4) make 

Finally create the executable:

`make`

## 5) run

Run the example models in the examples/ directory:

`./DRAGON examples/FILENAME.xml` 
m with installation!

## CREDITS

We acknowledge here the use of external routines/table:
* **dmspec.F** Routines for calculating the annihilation spectrum from [DarkSUSY](http://www.darksusy.org) package, to be cited as [Gondolo et al., 2004](http://arxiv.org/abs/astro-ph/0406204)
* **MilkyWay_DR0.5_DZ0.1_DPHI10_RMAX20_ZMAX5_galprop_format.fits.gz** The ISRF model used for the energy losses, to be downloaded from the [GALPROP](http://galprop.stanford.edu) package (v54) and to be cited as [Porter and Strong, 2008](http://adsabs.harvard.edu/abs/2008AAS...212.1810P)
* **webber_xsec.dat** Tabulated spallation cross sections, to be cited as [Webber et al., 2003](http://adsabs.harvard.edu/abs/2003ApJS..144..153W)
* **webber_xsec_total.dat** Tabulated inelastic cross sections, to be cited as [Webber et al., 2003](http://adsabs.harvard.edu/abs/2003ApJS..144..153W)
* **cparamlib** Library for secondary production in pp interactions from the [cparamlib](https://github.com/niklask/cparamlib) repository, to be cited as [Kamae, et al., 2007](https://arxiv.org/abs/astro-ph/0605581)
* **tinyxml** A C++ XML parser, from [here](http://www.grinninglizard.com/tinyxml)