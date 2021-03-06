#include "CSSHOWER++/Tools/Parton.H"
#include "CSSHOWER++/Tools/Singlet.H"
#include "ATOOLS/Phys/Cluster_Leg.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"

using namespace CSSHOWER;
using namespace ATOOLS;
using namespace std;

namespace CSSHOWER {
  std::ostream& operator<<(std::ostream& str, const Parton &part) {
    str<<"  Parton ["<<ATOOLS::ID(part.m_id)<<"], stat="
       <<part.m_stat<<", kin="<<part.m_kin<<", kscheme="<<part.m_kscheme
       <<", col="<<part.m_col<<" : "<<part.m_flav<<" : "<<part.m_mom
       <<" "<<sqrt(dabs(part.m_mom.Abs2()))<<" "<<sqrt(dabs(part.Mass2()))
       <<" ("<<part.GetFlow(1)<<","<<part.GetFlow(2)<<")"
       <<"["<<part.GetRFlow(1)<<","<<part.GetRFlow(2)<<"]"<<endl;
    if (part.m_pst==pst::IS)      str<<"     (Initial state parton)";
    else if (part.m_pst==pst::FS) str<<"     (Final state parton)  ";
    else                     str<<"                           ";
    str<<"  Colour partners ("
       <<(part.p_left?ATOOLS::ID(part.p_left->m_id):vector<int>())<<","
       <<(part.p_right?ATOOLS::ID(part.p_right->m_id):vector<int>())<<") "
       <<(part.p_colspec[0]?ATOOLS::ID(part.p_colspec[0]->m_id):vector<int>())<<","
       <<(part.p_colspec[1]?ATOOLS::ID(part.p_colspec[1]->m_id):vector<int>())<<endl;
    if (part.m_kt_soft[0]<std::numeric_limits<double>::max() ||
	part.m_kt_soft[1]<std::numeric_limits<double>::max()) {
      str<<"  k_T left : "<<sqrt(part.KtSoft(0))<<", k_T right : "<<sqrt(part.KtSoft(1))<<endl;
    }
    str<<"  k_T start : "<<sqrt(part.m_kt_start);
    str<<"  k_T test : "<<sqrt(part.m_kt_test);
    str<<"  k_T veto : "<<sqrt(part.m_kt_veto)<<"("<<sqrt(part.m_kt_max)<<")";
    str<<"  x_B : "<<part.m_xBj;
    str<<"  fromdec : "<<part.m_fromdec <<std::endl;
    if (part.p_prev || part.p_next) {
      if (part.p_prev) str<<"  P="<<part.p_prev;
      if (part.p_next) str<<"  N="<<part.p_next;
      str<<std::endl;
    }
    if (part.m_fixspec!=Vec4D())
      str<<"  fix spec : "<<part.m_fixspec<<", oldp : "<<part.OldMomentum()
	 <<" "<<(IsEqual(part.OldMomentum(),part.Momentum(),1.0e-6)?
		 "ok":"*** ERROR ***")<<"\n";
    return str;
  }
}

void Parton::DeleteAll()
{
  if (p_next) p_next->DeleteAll();
  delete this;
}

Parton *Parton::FollowUp()
{
  if (p_next) return p_next->FollowUp();
  return this;
}

bool Parton::Splits()
{
  if (this==p_sing->GetSplit()) return true;
  if (p_next == NULL) {
    return false;
  } else {
    return p_next->Splits();
  }
}

void Parton::UpdateDaughters()
{
  if (p_next==NULL) return;
  msg_Indent();
  msg_IODebugging()<<METHOD<<"("<<this<<") {\n";
  p_next->SetMomentum(m_mom);
  p_next->SetFlavour(m_flav);
  p_next->SetMass2(m_t);
  msg_IODebugging()<<*p_next;
  p_next->UpdateDaughters();
  msg_IODebugging()<<"}\n";
}

void Parton::UpdateNewDaughters(Parton *ref)
{
  if (p_next==NULL) return;
  ref=ref->GetSing()->GetLeft();
  if (ref==NULL) THROW(fatal_error,"Internal error");
  msg_Indent();
  msg_IODebugging()<<METHOD<<"("<<this<<") {\n";
  p_next->SetMomentum(m_mom);
  p_next->SetFlavour(m_flav);
  p_next->SetMass2(m_t);
  for (int n(1);n<=2;++n) {
    p_next->SetFlow(n,GetFlow(n));
    p_next->SetMEFlow(n,GetMEFlow(n));
  }
  p_next->SetStart(m_kt_start);
  p_next->SetKtMax(m_kt_max);
  p_next->SetVeto(m_kt_veto);
  p_next->SetId(m_id);
  p_next->SetFromDec(m_fromdec);
  msg_IODebugging()<<*p_next;
  p_next->UpdateNewDaughters(ref);
  msg_IODebugging()<<"}\n";
}

void Parton::UpdateColours()
{
  msg_IODebugging()<<METHOD<<"("<<this<<"): ("
		   <<GetMEFlow(1)<<","<<GetMEFlow(2)<<") -> ("
		   <<GetFlow(1)<<","<<GetFlow(2)<<") {\n";
  {
    msg_Indent();
    if (p_sing==NULL) THROW(fatal_error,"Cannot update flow");
    p_left=p_right=NULL;
    int f1(GetFlow(1)), f2(GetFlow(2));
    for (PLiter pit(p_sing->begin());pit!=p_sing->end();++pit) {
      if (f1 && f1==(*pit)->GetFlow(2)) (p_left=*pit)->SetRight(this);
      if (f2 && f2==(*pit)->GetFlow(1)) (p_right=*pit)->SetLeft(this);
    }
    msg_IODebugging()<<*this;
    if (this==p_sing->GetSplit() ||
	(p_prev && p_prev==p_sing->GetSplit())) {
      Parton *ds[2]={p_sing->GetLeft(),p_sing->GetRight()};
      for (int i(0);i<2;++i)
	for (int n(1);n<=2;++n)
	  if (ds[i]->GetFlow(n)==GetMEFlow(n)) {
	    ds[i]->SetFlow(n,GetFlow(n));
	    ds[i]->UpdateColours();
	  }
    }
    else if (p_next) {
      for (int n(1);n<=2;++n) p_next->SetFlow(n,GetFlow(n));
      p_next->UpdateColours();
    }
    for (int n(1);n<=2;++n) SetMEFlow(n,GetFlow(n));
  }
  msg_IODebugging()<<"}\n";
}

double Parton::Weight(const double &scale)
{
  double weight=1.0;
  for (size_t i(0);i<m_weights.size();++i)
    if (m_weights[i].first>scale) weight*=m_weights[i].second;
    else break;
  return weight;
}

void Parton::SetLeftOf(Parton * part)
{
  part->SetLeft(p_left);
  if (p_left) p_left->SetRight(part);
}

void Parton::SetRightOf(Parton * part)
{
  part->SetRight(p_right);
  if (p_right) p_right->SetLeft(part);
}

