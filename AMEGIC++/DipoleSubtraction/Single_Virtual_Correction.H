#ifndef AMEGIC_DipoleSubtraction_Single_Virtual_Correction_H
#define AMEGIC_DipoleSubtraction_Single_Virtual_Correction_H

#include "AMEGIC++/Main/Process_Base.H"

#include "AMEGIC++/DipoleSubtraction/Single_LOProcess.H"
#include "PHASIC++/Process/KP_Terms.H"

namespace PDF { class PDF_Base; } 

namespace PHASIC { class Virtual_ME2_Base; }

namespace AMEGIC {
  class Helicity;
  class Amplitude_Handler;
  class Phase_Space_Generator;
  class DipoleSplitting_Base;
  class Flavour_KernelsA;
  class Massive_KernelsA;

  class Single_Virtual_Correction : public Process_Base {
  protected:
    std::string             m_ptypename,m_pslibname;
    double                  m_iresult;
    double                  m_dalpha,m_dalpha_ff,
                            m_dalpha_fi,m_dalpha_if,m_dalpha_ii;
    Phase_Space_Generator * p_psgen;

    std::vector<std::vector<double> > p_dsij;

    std::vector<double>     m_cmur,m_wass;
    double                  m_x0,m_x1;

    Single_Virtual_Correction        * p_partner;

    Single_LOProcess      * p_LO_process;
    DipoleSplitting_Base  * p_dipole;
    PHASIC::Flavour_Kernels * p_flkern;
    PHASIC::Massive_Kernels * p_masskern;
    PHASIC::KP_Terms      * p_kpterms;
    PHASIC::Virtual_ME2_Base* p_loopme;  

    bool m_force_init, m_sccmur, m_kpcemode;
    bool m_checkborn, m_checkpoles, m_checkfinite;
    double m_checkpolesthreshold;
    int m_bvimode, m_user_bvimode, m_imode, m_checkloopmap;

    double m_bsum, m_vsum, m_isum, m_n;
    double m_mbsum, m_mvsum, m_misum, m_mn;
    double m_lastb, m_lastv, m_lasti, m_lastkp, m_lastki;


    std::vector<int> m_nomapflavs;
    bool m_loopmapped;
    /*------------------------------------------------------------------------------

      Constructors

      ------------------------------------------------------------------------------*/
  public:

    Single_Virtual_Correction();
    ~Single_Virtual_Correction();

    void AddPoint(const double &value);
    bool ReadIn(const std::string &pid);
    void WriteOut(const std::string &pid);
    void MPICollect(std::vector<double> &sv,size_t &i);
    void MPIReturn(std::vector<double> &sv,size_t &i);
    void MPISync(const int mode=0);

    bool FillIntegrator(PHASIC::Phase_Space_Handler *const psh);
    void SetScale(const PHASIC::Scale_Setter_Arguments &args);

    void SetSelector(const PHASIC::Selector_Key &key);
    void SetShower(PDF::Shower_Base *const ps);
    void SetFixedScale(const std::vector<double> &s);
    void SetSelectorOn(const bool on);
    void SetGenerator(PHASIC::ME_Generator_Base *const gen);
    ATOOLS::Flavour ReMap(const ATOOLS::Flavour &fl,const size_t &id) const;
    /*------------------------------------------------------------------------------

      Generic stuff for initialization of Single_Virtual_Correctiones

      ------------------------------------------------------------------------------*/
  protected:
    void                PolarizationNorm();
    double              Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
    bool                CreateChannelLibrary();

    /*------------------------------------------------------------------------------

      Initializing libraries, amplitudes, etc.

      ------------------------------------------------------------------------------*/
  public:
    void                AddChannels(std::list<std::string>*);
    bool                NewLibs();
    int                 InitAmplitude(Amegic_Model *,Topology *,
				      std::vector<Process_Base *> &,
				      std::vector<Process_Base *> &);
    bool                SetUpIntegrator();
    Single_LOProcess*   GetLOProcess()                 { return p_LO_process; }
    Amplitude_Handler * GetAmplitudeHandler()          { return p_partner->GetLOProcess()->GetAmplitudeHandler();}
    Helicity *          GetHelicity()                  { return p_partner->GetLOProcess()->GetHelicity(); }    
    double              Result()                       { return m_iresult; } 

    void                SelectLoopProcess();
    /*------------------------------------------------------------------------------

      Process management

      ------------------------------------------------------------------------------*/
  public:
    void             SetLookUp(const bool lookup);
    std::string      LibName()                          { return p_partner->GetLOProcess()->LibName();     }
    std::string      PSLibName()                        { return p_partner->GetLOProcess()->PSLibName();   }
    Process_Base   * Partner()       const              { return p_partner;     }
    void             Minimize();

    int Type() { return 200; }
    
    /*------------------------------------------------------------------------------

      Calculating total cross sections

      ------------------------------------------------------------------------------*/
  public:
    double         Partonic(const ATOOLS::Vec4D_Vector &,const int mode);
    double         DSigma(const ATOOLS::Vec4D_Vector &,bool,const int);
    double         operator()(const ATOOLS::Vec4D_Vector &,const int mode);
    void FillAmplitudes(std::vector<METOOLS::Spin_Amplitudes>& amps,
                        std::vector<std::vector<Complex> >& cols);

    int            NumberOfDiagrams();
    Point        * Diagram(int i);

    int PerformTests() { return 1; }

    bool Combinable(const size_t &idi,const size_t &idj);
    const ATOOLS::Flavour_Vector &CombinedFlavour(const size_t &idij);
    void FillProcessMap(PHASIC::NLOTypeStringProcessMap_Map *apmap);

    inline const PHASIC::KP_Terms*         KPTerms() { return p_kpterms; }
    inline const DipoleSplitting_Base*     Dipole()  { return p_dipole; }

  protected:
    double         Calc_Imassive(const ATOOLS::Vec4D *mom);
    double         Calc_I(const ATOOLS::Vec4D*);
    double         Calc_V_WhenMapped(const ATOOLS::Vec4D_Vector&);
    void           Calc_KP(const ATOOLS::Vec4D*, double, double, double, double, double);
    double         KPTerms(int mode, double scalefac2=1.0);
    double         Get_KPterms(PDF::PDF_Base *pdfa, PDF::PDF_Base *pdfb,
                               const double &, const double &,
                               ATOOLS::Flavour_Vector &,
                               double scalefac2=1.0);
    void           CheckPoleCancelation(const ATOOLS::Vec4D*);
    void           FillMEwgts(ATOOLS::ME_Weight_Info&);
    void           EndOptimize();

  };
}



#endif

