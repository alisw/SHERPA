pkglib_LTLIBRARIES       =	libSherpaMain.la

GITTAG = SHERPA/Main
include ../../git.make

libSherpaMain_la_SOURCES      =	Git_Info.C Sherpa.C Filter.C

localinc_HEADERS          =    Sherpa.H Filter.H

SYSLIBS	      =	$(FCLIBS) -lm -ldl $(CONDITIONAL_THREADLIBS)

libSherpaMain_la_LIBADD = -L$(SHERPABUILDDIR)/Initialization -L$(SHERPABUILDDIR)/Single_Events -L$(SHERPABUILDDIR)/Tools \
                -L$(SHERPABUILDDIR)/PerturbativePhysics -L$(SHERPABUILDDIR)/SoftPhysics -L$(SHERPABUILDDIR)/LundTools \
                -lSherpaInitialization -lSherpaSingleEvents \
                -lSherpaPerturbativePhysics -lSherpaSoftPhysics -lLundTools -lSherpaTools \
		@AMISICLIBS@ @CSSLIBS@ @DIRELIBS@ @ATOOLSLIBS@ @AHADICLIBS@ @SHRIMPSLIBS@ \
		@BEAMLIBS@ @COMIXLIBS@ @AMEGICLIBS@ @MODELLIBS@ @PDFLIBS@ \
		@MCATNLOLIBS@ @PHASICLIBS@ \
		@METOOLSLIBS@ @CONDITIONAL_HEPMC3LIBS@ @CONDITIONAL_HEPMC2LIBS@ \
		@HADRONSLIBS@ @PHOTONSLIBS@ @EXTRAXSLIBS@ \
                $(SYSLIBS)
