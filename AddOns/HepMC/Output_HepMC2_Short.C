#include "AddOns/HepMC/Output_HepMC2_Short.H"
#include "HepMC/GenEvent.h"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"

#ifdef USING__HEPMC2__IOGENEVENT
#include "HepMC/IO_GenEvent.h"
#endif

#ifdef USING__HEPMC2__DEFS
#include "HepMC/HepMCDefs.h"
#ifdef HEPMC_HAS_CROSS_SECTION
#include "HepMC/GenCrossSection.h"
#endif
#endif

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Output_HepMC2_Short::Output_HepMC2_Short(const Output_Arguments &args) :
  Output_Base("HepMC2S")
{
  m_basename=args.m_outpath+"/"+args.m_outfile;
  m_ext=".hepmc";
  int precision       = args.p_reader->GetValue<int>("OUTPUT_PRECISION",12);
#ifdef USING__HEPMC2__IOGENEVENT
  p_iogenevent = new HepMC::IO_GenEvent(*m_outstream.stream());
#ifdef HEPMC_HAS_CROSS_SECTION
  p_iogenevent->precision(precision);
#endif
#else
  THROW(fatal_error,"HepMC::IO_GenEvent asked for, but HepMC version too old.");
#endif
#ifdef HEPMC_HAS_CROSS_SECTION
  p_xs=new HepMC::GenCrossSection();
#endif
  p_event=new HepMC::GenEvent();
#ifdef USING__MPI
  if (mpi->Size()>1) {
    m_basename+="_"+rpa->gen.Variable("RNG_SEED");
  }
#endif
  m_outstream.open((m_basename+m_ext).c_str());
  if (!m_outstream.stream()->good())
    THROW(fatal_error, "Could not open event file "+m_basename+m_ext+".");
  m_outstream.stream()->precision(precision);
}

Output_HepMC2_Short::~Output_HepMC2_Short()
{
#ifdef USING__HEPMC2__IOGENEVENT
  delete p_iogenevent;
#endif
#ifdef HEPMC_HAS_CROSS_SECTION
  delete p_xs;
#endif
  m_outstream.close();
  delete p_event;
}

void Output_HepMC2_Short::SetXS(const double& xs, const double& xserr)
{
#ifdef HEPMC_HAS_CROSS_SECTION
  p_xs->set_cross_section(xs, xserr);
#endif
}

void Output_HepMC2_Short::Output(Blob_List* blobs, const double weight)
{
#ifdef USING__HEPMC2__IOGENEVENT
  p_event->clear();
  m_hepmc2.Sherpa2ShortHepMC(blobs, *p_event, weight);
  std::vector<HepMC::GenEvent*> subevents(m_hepmc2.GenSubEventList());
#ifdef HEPMC_HAS_CROSS_SECTION
  p_event->set_cross_section(*p_xs);
  for (size_t i(0); i<subevents.size(); ++i) {
    subevents[i]->set_cross_section(*p_xs);
  }
#endif
  if (subevents.size()) {
    for (size_t i(0); i<subevents.size(); ++i) {
      p_iogenevent->write_event(subevents[i]);
    }
    m_hepmc2.DeleteGenSubEventList();
  } else {
    p_iogenevent->write_event(p_event);
  }
#endif
}

void Output_HepMC2_Short::ChangeFile()
{
#ifdef USING__HEPMC2__IOGENEVENT
  delete p_iogenevent;
  m_outstream.close();
  std::string newname(m_basename+m_ext);
  for (size_t i(0);FileExists(newname);
       newname=m_basename+"."+ToString(++i)+m_ext);
  m_outstream.open(newname.c_str());
  if (!m_outstream.stream()->good())
    THROW(fatal_error, "Could not open event file "+newname+".");
  p_iogenevent = new HepMC::IO_GenEvent(*m_outstream.stream());
#endif
}

DECLARE_GETTER(Output_HepMC2_Short,"HepMC_Short",
	       Output_Base,Output_Arguments);

Output_Base *ATOOLS::Getter<Output_Base,Output_Arguments,Output_HepMC2_Short>::
operator()(const Output_Arguments &args) const
{
  return new Output_HepMC2_Short(args);
}

void ATOOLS::Getter<Output_Base,Output_Arguments,Output_HepMC2_Short>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"HepMC short output";
}

