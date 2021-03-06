#ifndef HADRONS_Current_Library_VA_P_S_H
#define HADRONS_Current_Library_VA_P_S_H

#include "HADRONS++/Current_Library/Current_Base.H"

namespace HADRONS
{

namespace VA_P_S_FFs
{
class FormFactor_Base {
protected:
  const ATOOLS::Flavour_Vector& m_flavs;
  const std::vector<int>& p_i;
  double m_m0, m_m1;
  double m_fplus, m_f0;
  bool   m_calced;
public:
  FormFactor_Base(struct GeneralModel model, double* masses, 
                  const ATOOLS::Flavour_Vector& flavs,
                  const std::vector<int>& indices) :
    m_flavs(flavs), p_i(indices){
    m_m0=masses[0];
    m_m1=masses[1];
    m_fplus=0.0; m_f0=0.0;
    m_calced=false;
  }
  virtual ~FormFactor_Base();
  virtual void CalcFFs( ATOOLS::Vec4D p0, ATOOLS::Vec4D p1 ) = 0;
  double fplus()  { if(!m_calced) abort(); return m_fplus; }
  double f0()     { if(!m_calced) abort(); return m_f0; }
};
}

class VA_P_S : public Current_Base {
  // 0=Scalar, 1=Pseudoscalar

  double m_Vxx;
  VA_P_S_FFs::FormFactor_Base* p_ff;
public:
  VA_P_S(const ATOOLS::Flavour_Vector& flavs,const std::vector<int>& indices,
         const std::string& name) :
    Current_Base(flavs, indices, name) {};
  ~VA_P_S() {
    if (p_ff!=NULL) delete p_ff;
  }
  void SetModelParameters( struct GeneralModel _md );
  void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);
};

}
#endif
