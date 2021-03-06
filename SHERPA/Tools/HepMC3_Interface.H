#ifndef SHERPA_Tools_HepMC3_Interface_H
#define SHERPA_Tools_HepMC3_Interface_H
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"

#ifdef USING__HEPMC3

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Phys/Weight_Info.H"
#include "PHASIC++/Process/Subprocess_Info.H"

#include "SHERPA/Tools/Variations.H"

#include <map>
#include <set>
#include <memory>

namespace HepMC3 {
  class GenEvent;
  class GenVertex;
  class GenRunInfo;
  class GenParticle;
  class GenCrossSection;
  class Writer;
  class WriterAscii;
  class WriterAsciiHepMC2;
  class WriterHEPEVT;
#ifdef USING__HEPMC3__WRITERROOTTREE
  class WriterRootTree;
#endif
#ifdef USING__HEPMC3__WRITERROOT
  class WriterRoot;
#endif

using GenVertexPtr = std::shared_ptr<GenVertex>;
using GenParticlePtr = std::shared_ptr<GenParticle>;

}

namespace HepMC = HepMC3;


namespace ATOOLS {
  class Blob_List;
  class Particle;
  class ME_Weight_Info;
  class NLO_subevtlist;
  class PDF_Info;
  class Variation_Weights;
}

namespace SHERPA {

  class EventInfo3 {
  private:
    ATOOLS::Blob * p_sp;
    bool m_usenamedweights, m_extendedweights;
    std::vector<SHERPA::Variations_Type::code> m_variationtypes;
    std::vector<double> m_orders;
    double m_wgt, m_mewgt, m_wgtnorm, m_ntrials, m_pswgt, m_pwgt;
    double m_mur2, m_muf12, m_muf22, m_muq2;
    double m_alphas, m_alpha;
    bool m_userhook;
    double m_userweight;
    PHASIC::nlo_type::code m_type;
    ATOOLS::ME_Weight_Info * p_wgtinfo;
    ATOOLS::PDF_Info * p_pdfinfo;
    ATOOLS::NLO_subevtlist * p_subevtlist;
    SHERPA::Variation_Weights * p_variationweights;

    void ReadIn(ATOOLS::Blob_Data_Base * &db, std::string name, bool abort);

  public:
    EventInfo3(ATOOLS::Blob * sp, const double &wgt,
              bool namedweights,
              bool extendedweights,
              bool includemeonlyweights);
    EventInfo3(const EventInfo3 &evtinfo);
    ~EventInfo3() {}
    bool WriteTo(HepMC::GenEvent &evt, const int& idx=-1);
    inline const ATOOLS::NLO_subevtlist * SubEvtList() { return p_subevtlist; }
    inline void SetWeight(const double &wgt)          { m_wgt=wgt; }
    inline void SetMEWeight(const double &mewgt)      { m_mewgt=mewgt; }
    inline void SetPartonicWeight(const double& pwgt) { m_pwgt=pwgt; }
    inline void SetMuR2(const double &mur2)           { m_mur2=mur2; }
    inline void SetMuF12(const double &muf12)         { m_muf12=muf12; }
    inline void SetMuF22(const double &muf22)         { m_muf22=muf22; }
    void SetAlphaS();
    void SetAlpha();
  };

  class HepMC3_Interface {
  private :
    bool m_usenamedweights, m_extendedweights, m_includemeonlyweights;
    bool m_hepmctree;

    std::map<ATOOLS::Blob*,HepMC::GenVertexPtr > m_blob2genvertex;
    std::map<ATOOLS::Particle*,HepMC::GenParticlePtr> m_particle2genparticle;
    std::set<ATOOLS::btp::code> m_ignoreblobs;

    HepMC::GenEvent *p_event;
    std::vector<HepMC::GenEvent*> m_subeventlist;

    bool Sherpa2ShortHepMC(const ATOOLS::Vec4D&,
                           const ATOOLS::Flavour&,
                           bool incoming,
                           HepMC::GenParticle*&);

  public :

    // constructor
    HepMC3_Interface();

    // destructor
    ~HepMC3_Interface();

    // member functions
    bool Sherpa2HepMC(ATOOLS::Blob_List *const blobs, HepMC::GenEvent& event);
    bool Sherpa2HepMC(ATOOLS::Blob_List *const blobs,                       std::shared_ptr<HepMC::GenRunInfo> run= std::shared_ptr<HepMC::GenRunInfo>());
    bool Sherpa2HepMCBlobtoGenVertex(ATOOLS::Blob *,     HepMC::GenVertexPtr &, HepMC::GenEvent& event ); 
    bool Sherpa2HepMC(ATOOLS::Particle *, HepMC::GenParticlePtr &);

    bool Sherpa2ShortHepMC(ATOOLS::Blob_List *const blobs, HepMC::GenEvent& event);
    bool Sherpa2ShortHepMC(ATOOLS::Blob_List *const blobs);
    bool SubEvtList2ShortHepMC(EventInfo3 &evtinfo);

    void AddCrossSection(HepMC::GenEvent& event,
                         const double& xs, const double& err);

    // inline functions
    inline HepMC::GenEvent* GenEvent() { return p_event; }
    // For RivetInterface to enable ME/shower only analysis of fully hadronised event
    inline void Ignore(ATOOLS::btp::code type) { m_ignoreblobs.insert(type); }

    inline HepMC::GenEvent* GenSubEvent(size_t i)
    { if (i<m_subeventlist.size()) return m_subeventlist[i]; return NULL; }
    inline std::vector<HepMC::GenEvent*> GenSubEventList()
    { return m_subeventlist; }

    inline void SetHepMCNamedWeights(bool x)   { m_usenamedweights=x; }
    inline bool HepMCNamedWeights()            { return m_usenamedweights; }
    inline void SetHepMCExtendedWeights(bool x){ m_extendedweights=x; }
    inline bool HepMCExtendedWeights()         { return m_extendedweights; }
    inline void SetHepMCTreeLike(bool x)       { m_hepmctree=x; }
    inline bool HepMCTreeLike()                { return m_hepmctree; }
    inline void SetHepMCIncludeMEOnlyVariations(bool x)
    { m_includemeonlyweights=x; }

    void DeleteGenSubEventList();

  };// end of class HepMC3_Interface

}// end of namespace SHERPA

#endif
#endif
