SUB_CPPFLAGS = \
	-I@abs_top_builddir@/include \
	-I@abs_top_builddir@/cparamlib

SUB_LINKFLAGS = \
    -L. -lDRAGON \
    -Lcparamlib -lcparamlib \
    -L. -lTiXML

LIBTOOL_DEPS = @LIBTOOL_DEPS@
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status --recheck

AM_CPPFLAGS = $(SUB_CPPFLAGS)
LDADD  = $(SUB_LNKFLAGS) 

SUBDIRS = cparamlib

export SUB_CPPFLAGS
export SUB_LNKFLAGS 
lib_LTLIBRARIES = libDRAGON.la libTiXML.la

libDRAGON_la_SOURCES = dragon.cc \
		 galaxy.cc \
         geometry.cc \
		 xsec.cc \
		 particle.cc \
		 crevolutor.cc \
		 eloss.cc \
         bfield.cc \
         diffusion.cc \
		 gas.cc \
		 grid.cc \
		 nucleilist.cc \
		 utilities.cc \
         kamae.cc \
         sources.cc \
         spectrum.cc \
		 input.cc \
		 dmspec.F \
		 WNEWTR_FUNC_aws.f \
		 YIELDX_011000_imos.f \
		 crn6.f \
		 include/constants.h \
         include/geometry.h \
         include/errorcode.h \
		 include/crevolutor.h \
		 include/galaxy.h \
		 include/diffusion.h \
		 include/gas.h \
		 include/grid.h \
		 include/input.h \
		 include/nucleilist.h \
		 include/particle.h \
		 include/utilities.h \
		 include/xsec.h \
		 include/kamae.h \
		 include/sources.h \
		 include/spectrum.h \
		 include/eloss.h \
         include/bfield.h \
		 include/dragon.h \
                 cparamlib/libcparamlib.la

libTiXML_la_SOURCES = tinystr.cc \
		 tinyxml.cc \
		 tinyxmlerror.cc \
		 tinyxmlparser.cc \
		 include/tinyxml.h \
		 include/tinystr.h 

bin_PROGRAMS = DRAGON
DRAGON_LDADD = libDRAGON.la libTiXML.la cparamlib/libcparamlib.la
$(SUB_LNKFLAGS)

DRAGON_SOURCES = main.cc \
    libDRAGON.la \
    libTiXML.la 

EXTRA_DIST = config_files/template.source.param \
	start.sh \
	cleancode.sh \
	output \
        m4 \
	data/isotope_list.dat \
	data/MilkyWay_DR0.5_DZ0.1_DPHI10_RMAX20_ZMAX5_galprop_format.fits.gz \
	data/galprop_barpol.dat \
	data/crxsecs_fragmentation_Evoli2019_cumulative_modified.dat \
	data/Winkler-ap_HeHe.dat \
	data/Winkler-ap_Hep.dat \
	data/Winkler-ap_pHe.dat \
	data/Winkler-ap_pp.dat \
	data/webber_xsec.dat \
	data/webber_xsec_total.dat \
	data/BUBBLE.DAT \
	examples/run_2D_Fluka_d033_1_1.xml \
	examples/run_2D_Fluka_d042_2.xml \
	examples/run_2D_Fluka_d050_1_1.xml \
	examples/run_2D_Fluka_d033_1_1.source.param \
	examples/run_2D_Fluka_d042_2.source.param \
	examples/run_2D_Fluka_d050_1_1.source.param \
	plots/plot_bc.pro \
	plots/plot_leptons.pro \
        plots/plot_leptons_DM.pro \
	plots/plot_protons.pro \
	plots/HEAO3_BC.txt \
	plots/PAMELA_posratio.txt \
	plots/PAMELA_protons.txt \
	plots/pyBC.py \
	plots/pyleptons.py \
	plots/pyprotons.py
	
dist-hook:
	rm -rf $(distdir)/output/* \
	rm -rf $(distdir)/*.la \
	rm -rf $(distdir)/m4/* 
