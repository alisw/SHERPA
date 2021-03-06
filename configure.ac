dnl Process this file with autoconf to produce a configure script.

AC_INIT(SHERPA MC, 2.2.8, sherpa@projects.hepforge.org, SHERPA-MC)
AC_CONFIG_MACRO_DIR([m4])
AM_CONFIG_HEADERS(config.h ATOOLS/Org/CXXFLAGS.H ATOOLS/Org/CXXFLAGS_PACKAGES.H)
AM_INIT_AUTOMAKE(tar-ustar)
AM_DISABLE_STATIC
AC_PREFIX_DEFAULT(${PWD})
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CANONICAL_BUILD
SHERPA_SETUP_BUILDSYSTEM
AC_DEFINE([SHERPA_NAME], ["Cho Oyu"], [Sherpa version name])
dnl Checks for programs.
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_CC

AC_PROG_FC(gfortran g95 g77 ifort)
AC_FC_LIBRARY_LDFLAGS
AC_FC_WRAPPERS
## Hack to make libtool use the FC Fortran compiler via F77, since
## libtool < 2.x doesn't understand --tag=FC
F77=$FC
AM_FFLAGS=$AM_FCFLAGS
AC_SUBST(F77)
AC_SUBST(AM_FFLAGS)

AC_ARG_ENABLE(mpi,
  AC_HELP_STRING([--enable-mpi], [Enable MPI]),
  [ AC_MSG_CHECKING(for MPI)
    case "${enableval}" in
      no)  AC_MSG_RESULT(no); mpi=false ;;
      yes) AC_MSG_RESULT(yes); mpi=true ;;
    esac ],
  [ AC_MSG_CHECKING(for MPI); AC_MSG_RESULT(no); mpi=false ])
if test "$mpi" = "true"; then
  AC_DEFINE([USING__MPI], "1", [using MPI])
  AX_PROG_CXX_MPI([test "$mpi" = "true"], [use_mpi=yes],
    [AC_MSG_FAILURE([MPI compiler not found.])])
  AS_AC_EXPAND(LOADMPIFORPY,from mpi4py import MPI)
else
  AC_PROG_CXX
  AS_AC_EXPAND(LOADMPIFORPY,)
fi
AM_CONDITIONAL(USING__MPI, test "$mpi" = "true" )

AC_ARG_ENABLE(pyext,
  AC_HELP_STRING([--enable-pyext], [Enable Python API]),
  [ AC_MSG_CHECKING(for Python extension)
    case "${enableval}" in
      no) AC_MSG_RESULT(no); pyext=false ;;
      yes)  AC_MSG_RESULT(yes); pyext=true ;;
    esac ],
  [ AC_MSG_CHECKING(for Python extension); AC_MSG_RESULT(no); pyext=false])
if test x$pyext == xtrue; then
  AM_PATH_PYTHON
  AX_PYTHON_DEVEL
  AX_PKG_SWIG(1.3.17,[],[ AC_MSG_ERROR([SWIG is required to build..]) ])
  AX_SWIG_ENABLE_CXX
  AX_SWIG_MULTI_MODULE_SUPPORT
  AX_SWIG_PYTHON
fi
AM_CONDITIONAL(ENABLE_PYEXT, [test x$pyext == xtrue])

AC_ARG_ENABLE(ufo,
  AC_HELP_STRING([--enable-ufo], [Enable UFO support]),
  [ AC_MSG_CHECKING(for UFO support)
    case "${enableval}" in
      no) AC_MSG_RESULT(no); ufo=false ;;
      yes)  AC_MSG_RESULT(yes); ufo=true ;;
    esac ],
  [ AC_MSG_CHECKING(for UFO support); AC_MSG_RESULT(no); ufo=false])
if test x$ufo == xtrue; then
  AM_PATH_PYTHON
fi
AM_CONDITIONAL(ENABLE_UFO, [test x$ufo == xtrue])

case "$build_os:$build_cpu:$build_vendor" in
  *darwin*:*:*) ;;
  *linux*:*:*) CXX="$CXX -Wl,--no-as-needed" ;;
  *) CXX="$CXX -Wl,--no-as-needed" ;;
esac

AM_PROG_LIBTOOL
dnl Checks for header files.
AC_HEADER_STDC
dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
dnl Checks for library functions.
AC_CHECK_FUNCS(mkdir)

AX_LIB_SQLITE3(3.3.9)
AS_IF([test "x$SQLITE3_VERSION" = x],
  [ AC_MSG_ERROR([Sqlite3 required but not found.
To install a local copy, use the configure option --with-sqlite3=install.]) ])

SHERPA_SETUP_CONFIGURE_OPTIONS
SHERPA_SETUP_VARIABLES

AC_OUTPUT(ATOOLS/Math/Makefile \
	  ATOOLS/Phys/Makefile \
	  ATOOLS/Org/Makefile \
	  ATOOLS/Makefile \
          METOOLS/Explicit/Makefile \
          METOOLS/Currents/Makefile \
          METOOLS/Vertices/Makefile \
          METOOLS/Colors/Makefile \
          METOOLS/SpinCorrelations/Makefile \
          METOOLS/Main/Makefile \
          METOOLS/Loops/Makefile \
          METOOLS/Makefile \
	  BEAM/Main/Makefile \
	  BEAM/Makefile \
	  PDF/Main/Makefile \
	  PDF/Electron/Makefile \
	  PDF/LHAPDF/Makefile \
          PDF/CT12/Makefile \
          PDF/CT10/Makefile \
          PDF/MSTW/Makefile \
	  PDF/MRST/Makefile \
	  PDF/NNPDF/Makefile \
	  PDF/CTEQ/Makefile \
	  PDF/GRV/Makefile \
	  PDF/Remnant/Makefile \
	  PDF/Makefile \
	  MODEL/Main/Makefile \
	  MODEL/SM/Makefile \
	  MODEL/HEFT/Makefile \
	  MODEL/TauPi/Makefile \
	  MODEL/UFO/Makefile \
	  MODEL/UFO/Sherpa-generate-model \
	  MODEL/UFO/sconstruct_template \
	  MODEL/Makefile \
	  MCATNLO/Tools/Makefile \
	  MCATNLO/Calculators/Makefile \
	  MCATNLO/Showers/Makefile \
	  MCATNLO/Main/Makefile \
	  MCATNLO/Makefile \
	  PHASIC++/Main/Makefile \
	  PHASIC++/Channels/Makefile \
	  PHASIC++/Process/Makefile \
	  PHASIC++/Selectors/Makefile \
	  PHASIC++/Scales/Makefile \
	  PHASIC++/Enhance/Makefile \
	  PHASIC++/Decays/Makefile \
	  PHASIC++/Makefile \
	  EXTRA_XS/Main/Makefile \
	  EXTRA_XS/Two2Two/Makefile \
	  EXTRA_XS/NLO/Makefile \
	  EXTRA_XS/Cluster/Makefile \
          EXTRA_XS/One2Two/Makefile \
	  EXTRA_XS/One2Three/Makefile \
	  EXTRA_XS/Makefile \
	  AMEGIC++/Main/Makefile \
          AMEGIC++/Main/makelibs \
          AMEGIC++/Main/make2scons \
	  AMEGIC++/DipoleSubtraction/Makefile \
          AMEGIC++/Amplitude/Makefile \
          AMEGIC++/Amplitude/Zfunctions/Makefile \
          AMEGIC++/Phasespace/Makefile \
          AMEGIC++/String/Makefile \
          AMEGIC++/Cluster/Makefile \
          AMEGIC++/Makefile \
	  COMIX/Amplitude/Makefile \
	  COMIX/Phasespace/Makefile \
	  COMIX/Main/Makefile \
	  COMIX/Cluster/Makefile \
	  COMIX/Makefile \
	  CSSHOWER++/Tools/Makefile \
	  CSSHOWER++/Calculators/Makefile \
	  CSSHOWER++/Showers/Makefile \
	  CSSHOWER++/Main/Makefile \
	  CSSHOWER++/Makefile \
	  DIRE/Tools/Makefile \
	  DIRE/Shower/Makefile \
	  DIRE/Gauge/Makefile \
	  DIRE/Lorentz/Makefile \
	  DIRE/Main/Makefile \
	  DIRE/Makefile \
	  AMISIC++/Tools/Makefile \
	  AMISIC++/Model/Makefile \
	  AMISIC++/Main/Makefile \
	  AMISIC++/Makefile \
	  AHADIC++/Tools/Makefile \
	  AHADIC++/Formation/Makefile \
	  AHADIC++/Decays/Makefile \
	  AHADIC++/Main/Makefile \
	  AHADIC++/Makefile \
	  HADRONS++/Main/Makefile \
	  HADRONS++/ME_Library/Makefile \
	  HADRONS++/Current_Library/Makefile \
	  HADRONS++/PS_Library/Makefile \
	  HADRONS++/Run/Makefile \
	  HADRONS++/Makefile \
	  PHOTONS++/Main/Makefile \
          PHOTONS++/Tools/Makefile \
          PHOTONS++/PhaseSpace/Makefile \
          PHOTONS++/MEs/Makefile \
	  PHOTONS++/Makefile \
	  SHRiMPS/Makefile \
	  SHRiMPS/Main/Makefile \
	  SHRiMPS/Event_Generation/Makefile \
	  SHRiMPS/Beam_Remnants/Makefile \
	  SHRiMPS/Cross_Sections/Makefile \
	  SHRiMPS/Eikonals/Makefile \
	  SHRiMPS/Tools/Makefile \
	  SHERPA/Tools/Makefile \
	  SHERPA/LundTools/Makefile \
	  SHERPA/LundTools/pythia-6.4.18.f \
	  SHERPA/Initialization/Makefile \
	  SHERPA/PerturbativePhysics/Makefile \
	  SHERPA/SoftPhysics/Makefile \
	  SHERPA/Single_Events/Makefile \
	  SHERPA/Main/Makefile \
	  SHERPA/Run/Makefile \
          SHERPA/Run/Sherpa-config \
	  SHERPA/Makefile \
	  Manual/Makefile \
          Manual/sherpa-completion \
	  Examples/Makefile \
	  Examples/API/ME2-Python/test.py \
	  Examples/API/Events/test.py \
	  Examples/API/MPIEvents/test.py \
	  AddOns/Python/Makefile \
	  AddOns/Analysis/Tools/Makefile \
	  AddOns/Analysis/Main/Makefile \
	  AddOns/Analysis/Triggers/Makefile \
	  AddOns/Analysis/Observables/Makefile \
	  AddOns/Analysis/Analyses/Makefile \
	  AddOns/Analysis/Scripts/Makefile \
	  AddOns/Analysis/Makefile \
	  AddOns/Pythia/Makefile \
	  AddOns/Higgs/Makefile \
	  AddOns/BlackHat/Makefile \
	  AddOns/OpenLoops/Makefile \
	  AddOns/Recola/Makefile \
          AddOns/MCFM/Makefile \
          AddOns/LH_OLE/Makefile \
	  AddOns/HZTool/Makefile \
	  AddOns/PGS/Makefile \
	  AddOns/Delphes/Makefile \
	  AddOns/Rivet/Makefile \
	  AddOns/HepMC/Makefile \
	  AddOns/Root/Makefile \
	  AddOns/Makefile \
	  Makefile )

if test "x$VERSIONING" != "x"; then
  mv "SHERPA/Run/Sherpa-config" "SHERPA/Run/Sherpa$VERSIONING-config"
fi
