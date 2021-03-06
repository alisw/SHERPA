#ifndef CSSHOWER_Showers_Shower_H
#define CSSHOWER_Showers_Shower_H

#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "CSSHOWER++/Showers/Sudakov.H"
#include "CSSHOWER++/Tools/Singlet.H"
#include "CSSHOWER++/Showers/Kinematics_Base.H"

namespace SHERPA { class Variation_Weights; }

namespace ATOOLS { class Cluster_Leg; }

namespace CSSHOWER {

  class Shower {
  private:
    Kinematics_FF      m_kinFF;
    Kinematics_FI      m_kinFI;
    Kinematics_IF      m_kinIF;
    Kinematics_II      m_kinII;
    ATOOLS::Flavour    m_flavA, m_flavB, m_flavC;
    Singlet      *     p_actual;
    Sudakov            m_sudakov;
    PDF::ISR_Handler * p_isr;
    Parton           * m_last[4];
    int                m_kscheme, m_noem, m_recdec;
    bool               m_use_bbw;

    ATOOLS::Cluster_Leg *p_old[2];

    std::map<std::string,double,ATOOLS::String_Sort> m_efac;

    double m_weight;
    bool m_reweight;
    SHERPA::Variation_Weights * p_variationweights;
    double m_maxreweightfactor;

    void   ResetScales(const double &kt2);
    void   SetSplitInfo(const ATOOLS::Vec4D &psplit,
			const ATOOLS::Vec4D &pspect,
			Parton *const split,Parton *const newb,
			Parton *const newc,const int mode);
    double GetXBj(Parton *const p) const;
    int    SetXBj(Parton *const p) const;
    Parton *SelectSplitting(double &);
    bool   EvolveSinglet(Singlet *,const size_t &,size_t &);
    bool   PerformSplitting();
    bool   TrialEmission(double &,Parton *);
    int    RemnantTest(Parton *const p,const ATOOLS::Poincare_Sequence *lt);
    int    UpdateDaughters(Parton *const split,Parton *const newpB,
			   Parton *const newpC,int mode);
    int MakeKinematics(Parton *const split,const ATOOLS::Flavour &fla,
		       const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc,
		       const int mode,const int fc);

    ATOOLS::Cluster_Amplitude *  SingletToAmplitude(const Singlet * sing, ATOOLS::Cluster_Amplitude * ampl);
    void   PartonToAmplitude(const Parton *  parton, ATOOLS::Cluster_Amplitude *);
    void   CheckAmplitude(const ATOOLS::Cluster_Amplitude *ampl);
 public:
    Shower(PDF::ISR_Handler * isr,const int qed,
	   ATOOLS::Data_Reader *const dataread,int type);
    ~Shower();
    bool   EvolveShower(Singlet *,const size_t &,size_t &);
    int    ReconstructDaughters(Singlet *const split,const int mode,
				Parton *const pi=NULL,Parton *const pj=NULL);
    double Reweight(SHERPA::Variation_Parameters*,
                    SHERPA::Variation_Weights*,
                    Parton& splitter);
    void SetMS(ATOOLS::Mass_Selector *const ms);
    double EFac(const std::string &sfk) const;
    inline Sudakov *GetSudakov() { return &m_sudakov; }
    inline const Sudakov *GetSudakov() const { return &m_sudakov; }

    inline const Kinematics_FF *KinFF() const { return &m_kinFF; }
    inline const Kinematics_FI *KinFI() const { return &m_kinFI; }
    inline const Kinematics_IF *KinIF() const { return &m_kinIF; }
    inline const Kinematics_II *KinII() const { return &m_kinII; }

    inline Kinematics_FF *KinFF() { return &m_kinFF; }
    inline Kinematics_FI *KinFI() { return &m_kinFI; }
    inline Kinematics_IF *KinIF() { return &m_kinIF; }
    inline Kinematics_II *KinII() { return &m_kinII; }

    inline double Weight() const { return m_weight; }

    inline int KinScheme() const { return m_kscheme; }

    inline PDF::ISR_Handler *ISR() { return p_isr; }

    inline void SetVariationWeights(SHERPA::Variation_Weights * varwgts)
    { p_variationweights = varwgts; }

    inline void SetUsesBBW(bool on) { m_use_bbw = on; }

  };
}

#endif
