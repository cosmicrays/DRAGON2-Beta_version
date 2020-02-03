lib_LTLIBRARIES = libcparamlib.la
#DRAGON_LDADD = libDRAGON.la libTiXML.la libcparamlib.la

libcparamlib_la_SOURCES = \
    antinumu.c \
    elec.c \
    nue.c \
    posi.c \
    antinue.c \
    cparamlib.h \
    gamma_param.c \
    numu.c \
    sigma.c

#EXTRA_DIST = doc/Doxyfile
dist-hook:
	rm -rf `find $(distdir) -type d -name CVS`
