#include "PHASIC++/Process/Process_Base.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Channels/BBar_Multi_Channel.H"
#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Decay_Info.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/My_MPI.H"
#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <algorithm>

using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

int Process_Base::s_usefmm(-1);

Process_Base::Process_Base():
  p_parent(NULL), p_selected(this), p_mapproc(NULL),
  p_sproc(NULL), p_proc(this),
  p_int(new Process_Integrator(this)), p_selector(NULL),
  p_cuts(NULL), p_gen(NULL), p_shower(NULL), p_mc(NULL),
  p_scale(NULL), p_kfactor(NULL),
  m_nin(0), m_nout(0), m_maxcpl(2,99), m_mincpl(2,0), 
  m_tinfo(1), m_mcmode(0), m_cmode(0),
  m_lookup(false), m_use_biweight(true), p_apmap(NULL),
  p_variationweights(NULL), m_variationweightsowned(false),
  p_lkfvariationweights(NULL)
{
  m_last=m_lastb=0.0;
  if (s_usefmm<0) s_usefmm=ToType<int>(rpa->gen.Variable("PB_USE_FMM"));
}

Process_Base::~Process_Base() 
{
  if (p_kfactor) delete p_kfactor;
  if (p_scale) delete p_scale;
  if (m_variationweightsowned && p_variationweights) delete p_variationweights;
  if (p_lkfvariationweights) delete p_lkfvariationweights;
  delete p_selector;
  delete p_int;
}

Process_Base *Process_Base::Selected()
{ 
  if (!p_selected) return NULL;
  if (p_selected!=this) return p_selected->Selected();
  return this; 
}

bool Process_Base::SetSelected(Process_Base *const proc)
{
  if (proc==this) {
    p_selected=this;
    return true;
  }
  if (IsGroup())
    for (size_t i(0);i<Size();++i)
      if ((*this)[i]->SetSelected(proc)) {
	p_selected=(*this)[i];
	return true;
      }
  return false;
}

Process_Base *Process_Base::Parent()
{ 
  if (p_parent && p_parent!=this) return p_parent->Parent();
  return this; 
}

bool Process_Base::GeneratePoint()
{
  return true;
}

void Process_Base::AddPoint(const double &value)
{
}

bool Process_Base::ReadIn(const std::string &pid)
{
  return true;
}

void Process_Base::WriteOut(const std::string &pid)
{
}

void Process_Base::EndOptimize()
{
}

void Process_Base::MPICollect(std::vector<double> &sv,size_t &i)
{
  if (IsGroup())
    for (size_t j(0);j<Size();++j)
      (*this)[j]->MPICollect(sv,i);
}

void Process_Base::MPIReturn(std::vector<double> &sv,size_t &i)
{
  if (IsGroup())
    for (size_t j(0);j<Size();++j)
      (*this)[j]->MPIReturn(sv,i);
}

void Process_Base::MPISync(const int mode)
{
  if (mode) return;
#ifdef USING__MPI
  size_t i(0), j(0);
  std::vector<double> sv;
  MPICollect(sv,i);
  if (mpi->Size()>1)
    mpi->Allreduce(&sv[0],sv.size(),MPI_DOUBLE,MPI_SUM);
  MPIReturn(sv,j);
#endif
}

void Process_Base::SetFixedScale(const std::vector<double> &s)
{
  if (IsMapped()) p_mapproc->SetFixedScale(s);
  if (p_scale!=NULL) p_scale->SetFixedScale(s);
}

void Process_Base::SetSelectorOn(const bool on)
{
  Selector()->SetOn(on);
}

void Process_Base::SetUseBIWeight(bool on)
{
  m_use_biweight=on;
}

double Process_Base::Differential(const Cluster_Amplitude &ampl,int mode) 
{
  Vec4D_Vector p(ampl.Legs().size());
  for (size_t i(0);i<ampl.NIn();++i) p[i]=-ampl.Leg(i)->Mom();
  if (mode&16) THROW(not_implemented,"Invalid mode");
  for (size_t i(ampl.NIn());i<p.size();++i) p[i]=ampl.Leg(i)->Mom();
  if (mode&64) {
    if (mode&1) return 1.0;
    return Trigger(p);
  }
  bool selon(Selector()->On());
  if (!Trigger(p)) {
    if ((mode&1) && selon) {
      SetSelectorOn(false);
      Trigger(p);
    }
  }
  if (mode&2) {
    std::vector<double> s(ScaleSetter(1)->Scales().size(),0.0);
    s[stp::fac]=ampl.MuF2();
    s[stp::ren]=ampl.MuR2();
    s[stp::res]=ampl.MuQ2();
    if (s.size()>stp::size+stp::res)
      s[stp::size+stp::res]=ampl.KT2();
    SetFixedScale(s);
  }
  if (mode&4) SetUseBIWeight(false);
  if (mode&128) this->GeneratePoint(); 
  double res(this->Differential(p));
  if (mode&32) {
    SP(Phase_Space_Handler) psh(Parent()->Integrator()->PSHandler());
    res*=psh->Weight(p);
  }
  if (mode&4) SetUseBIWeight(true);
  if (mode&2) SetFixedScale(std::vector<double>());
  if (Selector()->On()!=selon) SetSelectorOn(selon);
  return res;
}

bool Process_Base::IsGroup() const
{
  return false;
}

bool Process_Base::FillIntegrator
(Phase_Space_Handler *const psh)
{
  return false;
}

bool Process_Base::InitIntegrator
(Phase_Space_Handler *const psh)
{
  if (!p_sproc) return true;
  DEBUG_FUNC(m_name);
  SetBBarMC(new BBar_Multi_Channel(this,p_sproc,psh));
  psh->SetFSRIntegrator(p_mc);
  return true;
}

void Process_Base::UpdateIntegrator
(Phase_Space_Handler *const psh)
{
}

class Order_Flavour {
  FMMap* p_fmm;
  int Order_SVFT(const Flavour &a,const Flavour &b) 
  {
    if (a.IsScalar() && !b.IsScalar()) return 1;
    if (a.IsVector() && !b.IsScalar() && 
	!b.IsVector()) return 1;
    if (a.IsFermion() && !b.IsFermion() && 
	!b.IsScalar() && !b.IsVector()) return 1;
    return 0;
  }
  int Order_Multi(const Flavour &a,const Flavour &b)
  {
    if ((*p_fmm)[int(a.Kfcode())]==0 || 
	(*p_fmm)[int(b.Kfcode())]==0) return 0;
    if ((*p_fmm)[int(a.Kfcode())]>
	(*p_fmm)[int(b.Kfcode())]) return 1;
    return 0;
  }
  int operator()(const Flavour &a,const Flavour &b)
  {
    if (a.Priority()>b.Priority()) return 1;
    if (a.Priority()<b.Priority()) return 0;
    if (!a.Strong()&&b.Strong()) return 1;
    if (a.Strong()&&!b.Strong()) return 0;
    if (p_fmm) {
      if (Order_Multi(a,b)) return 1;
      if (Order_Multi(b,a)) return 0;
    }
    if (a.Mass()>b.Mass()) return 1;
    if (a.Mass()<b.Mass()) return 0;
    if (Order_SVFT(a,b)) return 1;
    if (Order_SVFT(b,a)) return 0;
    if (!a.IsAnti()&&b.IsAnti()) return 1;
    if (a.IsAnti()&&!b.IsAnti()) return 0;
    return a.Kfcode()<b.Kfcode();
  }
public:
  Order_Flavour(FMMap* fmm): p_fmm(fmm) {}
  int operator()(const Subprocess_Info &a,const Subprocess_Info &b)
  { return (*this)(a.m_fl,b.m_fl); }
  int operator()(const Cluster_Leg *a,const Cluster_Leg *b)
  { return (*this)(a->Flav(),b->Flav()); }
};// end of class Order_Flavour

void Process_Base::SortFlavours(Subprocess_Info &info,FMMap *const fmm)
{
  if (info.m_ps.empty()) return;
  ATOOLS::Flavour heaviest(kf_photon);
  for (size_t i(0);i<info.m_ps.size();++i) {
    if (info.m_ps[i].m_fl.Mass()>heaviest.Mass()) heaviest=info.m_ps[i].m_fl;
    else if (info.m_ps[i].m_fl.Mass()==heaviest.Mass() &&
	     !info.m_ps[i].m_fl.IsAnti()) heaviest=info.m_ps[i].m_fl;
  }
  std::sort(info.m_ps.begin(),info.m_ps.end(),Order_Flavour(fmm));
  for (size_t i(0);i<info.m_ps.size();++i) SortFlavours(info.m_ps[i]);
}

void Process_Base::SortFlavours(Process_Info &pi,const int mode)
{
  FMMap fmm;
  for (size_t i(0);i<pi.m_ii.m_ps.size();++i) {
    const Flavour *hfl=&pi.m_ii.m_ps[i].m_fl;
    if (fmm.find(int(hfl->Kfcode()))==fmm.end()) 
      fmm[int(hfl->Kfcode())]=0;
    if (hfl->IsFermion()) {
      fmm[int(hfl->Kfcode())]+=10;
      if (!hfl->IsAnti()) fmm[int(hfl->Kfcode())]+=10;
    }
  }
  for (size_t i(0);i<pi.m_fi.m_ps.size();++i) {
    const Flavour *hfl=&pi.m_fi.m_ps[i].m_fl;
    if (fmm.find(int(hfl->Kfcode()))==fmm.end()) 
      fmm[int(hfl->Kfcode())]=0;
    if (hfl->IsFermion()) fmm[int(hfl->Kfcode())]++;
  }
  if (mode&1) SortFlavours(pi.m_ii,s_usefmm?&fmm:NULL);
  SortFlavours(pi.m_fi,s_usefmm?&fmm:NULL);
}

void Process_Base::Init(const Process_Info &pi,
			BEAM::Beam_Spectra_Handler *const beamhandler,
			PDF::ISR_Handler *const isrhandler,const int mode)
{
  m_pinfo=pi;
  m_nin=m_pinfo.m_ii.NExternal();
  m_nout=m_pinfo.m_fi.NExternal();
  m_flavs.resize(m_nin+m_nout);
  if (m_pinfo.m_ii.m_ps.size()>0 && m_pinfo.m_fi.m_ps.size()>0) {
    if (!(mode&1)) SortFlavours(m_pinfo);
    std::vector<Flavour> fl;
    m_pinfo.m_ii.GetExternal(fl);
    m_pinfo.m_fi.GetExternal(fl);
    if (fl.size()!=m_nin+m_nout) THROW(fatal_error,"Internal error");
    for (size_t i(0);i<fl.size();++i) m_flavs[i]=fl[i];
    m_name=GenerateName(m_pinfo.m_ii,m_pinfo.m_fi);
    m_pinfo.m_fi.BuildDecayInfos(m_nin);
    m_decins=m_pinfo.m_fi.GetDecayInfos();
    if (IsGroup()) {
      if (m_pinfo.m_nminq>0 || m_pinfo.m_nmaxq<m_nin+m_nout) 
        m_name+="__NQ_"+ToString(m_pinfo.m_nminq)+
	  "-"+ToString(m_pinfo.m_nmaxq);
    }
  }
  double massin=0.0, massout=0.0;
  for (size_t i=0;i<m_nin;++i) massin+=m_flavs[i].Mass();
  for (size_t i=m_nin;i<m_nin+m_nout;++i) massout+=m_flavs[i].Mass();
  p_int->SetISRThreshold(Max(massin,massout));
  p_int->Initialize(beamhandler,isrhandler);
  m_symfac=m_pinfo.m_fi.FSSymmetryFactor();
  m_name+=pi.m_addname;
}

std::string Process_Base::GenerateName(const Subprocess_Info &info) 
{
  std::string name(info.m_fl.IDName());
  if (info.m_fl.Kfcode()==kf_quark && info.m_fl.IsAnti()) name+="b";
  if (info.m_ps.empty()) return name;
  name+="["+GenerateName(info.m_ps.front());
  for (size_t i(1);i<info.m_ps.size();++i) 
    name+="__"+GenerateName(info.m_ps[i]);
  if (info.m_nloqcdtype!=nlo_type::lo) 
    name+="__QCD("+ToString(info.m_nloqcdtype)+info.m_sv+")";
  if (info.m_nloewtype!=nlo_type::lo) 
    name+="__EW("+ToString(info.m_nloewtype)+info.m_sv+")";
  return name+="]";
}

std::string Process_Base::GenerateName
(const Subprocess_Info &ii,const Subprocess_Info &fi) 
{
  char nii[3], nfi[3];
  if (sprintf(nii,"%i",(int)ii.NExternal())<=0)
    THROW(fatal_error,"Conversion error");
  if (sprintf(nfi,"%i",(int)fi.NExternal())<=0)
    THROW(fatal_error,"Conversion error");
  std::string name(nii+std::string("_")+nfi);
  for (size_t i(0);i<ii.m_ps.size();++i) name+="__"+GenerateName(ii.m_ps[i]);
  for (size_t i(0);i<fi.m_ps.size();++i) name+="__"+GenerateName(fi.m_ps[i]);
  if (fi.m_nloqcdtype!=nlo_type::lo) 
    name+="__QCD("+ToString(fi.m_nloqcdtype)+fi.m_sv+")";
  if (fi.m_nloewtype!=nlo_type::lo) 
    name+="__EW("+ToString(fi.m_nloewtype)+fi.m_sv+")";
  return name;
}

class Order_NDecay {
public:
  int operator()(const Decay_Info *a,const Decay_Info *b) 
  { return IdCount(a->m_id)>IdCount(b->m_id); }
};// end of class Order_NDecay

void Process_Base::SortFlavours
(std::vector<Cluster_Leg*> &legs,FMMap *const fmm)
{
  if (legs.empty()) return;
  ATOOLS::Flavour heaviest(kf_photon);
  for (size_t i(0);i<legs.size();++i) {
    if (legs[i]->Flav().Mass()>heaviest.Mass()) heaviest=legs[i]->Flav();
    else if (legs[i]->Flav().Mass()==heaviest.Mass() &&
	     !legs[i]->Flav().IsAnti()) heaviest=legs[i]->Flav();
  }
  std::sort(legs.begin(),legs.end(),Order_Flavour(fmm));
}

void Process_Base::SortFlavours
(Cluster_Amplitude *const ampl,const int mode)
{
  FMMap fmm;
  DecayInfo_Vector cs;
  ClusterLeg_Vector il, fl;
  std::vector<int> dec(ampl->Legs().size(),0);
  std::map<size_t,ClusterLeg_Vector> dmap;
  for (size_t j(0);j<ampl->Decays().size();++j) {
    Decay_Info *cdi(ampl->Decays()[j]);
    size_t did(cdi->m_id), ndc(IdCount(did));
    for (size_t i(ampl->NIn());i<dec.size();++i)
      if (did&ampl->Leg(i)->Id()) {
	dec[i]=1;
	dmap[cdi->m_id].push_back(ampl->Leg(i));
      }
    bool core(true);
    for (size_t i(0);i<ampl->Decays().size();++i)
      if ((ampl->Decays()[i]->m_id&did) &&
	  IdCount(ampl->Decays()[i]->m_id)>ndc) {
	core=false;
	break;
      }
    if (!core) continue;
    int kfc(cdi->m_fl.Kfcode());
    if (fmm.find(kfc)==fmm.end()) fmm[kfc]=0;
    if (cdi->m_fl.IsFermion()) {
      fmm[kfc]+=10;
      if (!cdi->m_fl.IsAnti()) fmm[kfc]+=10;
    }
    cs.push_back(cdi);
  }
  for (size_t i(0);i<ampl->Legs().size();++i)
    if (i<ampl->NIn()) {
      ampl->Leg(i)->SetFlav(ampl->Leg(i)->Flav().Bar());
      il.push_back(ampl->Leg(i));
      int kfc(ampl->Leg(i)->Flav().Kfcode());
      if (fmm.find(kfc)==fmm.end()) fmm[kfc]=0;
      if (ampl->Leg(i)->Flav().IsFermion()) {
	fmm[kfc]+=10;
	if (!ampl->Leg(i)->Flav().IsAnti()) fmm[kfc]+=10;
      }
    }
    else {
      if (dec[i]) continue;
      fl.push_back(ampl->Leg(i));
      int kfc(ampl->Leg(i)->Flav().Kfcode());
      if (fmm.find(kfc)==fmm.end()) fmm[kfc]=0;
      if (ampl->Leg(i)->Flav().IsFermion()) ++fmm[kfc];
    }
  if (mode&1) SortFlavours(il,s_usefmm?&fmm:NULL);
  for (size_t i(0);i<cs.size();++i) {
    ampl->CreateLeg(Vec4D(),cs[i]->m_fl,ColorID(),cs[i]->m_id);
    fl.push_back(ampl->Legs().back());
  }
  SortFlavours(fl,s_usefmm?&fmm:NULL);
  if (cs.size()) {
    cs=ampl->Decays();
    std::sort(cs.begin(),cs.end(),Order_NDecay());
    while (fl.size()<ampl->Legs().size()-ampl->NIn())
      for (ClusterLeg_Vector::iterator
	     fit(fl.begin());fit!=fl.end();++fit)
	if (dmap.find((*fit)->Id())!=dmap.end()) {
	  ClusterLeg_Vector cl(dmap[(*fit)->Id()]);
	  size_t inc(0), ncd(IdCount((*fit)->Id()));
	  for (size_t i(0);i<cs.size();++i)
	    if (IdCount(cs[i]->m_id)<ncd &&
		(cs[i]->m_id&(*fit)->Id()) && (cs[i]->m_id&inc)==0) {
	    ampl->CreateLeg(Vec4D(),cs[i]->m_fl,ColorID(),cs[i]->m_id);
	    for (ClusterLeg_Vector::iterator
		   cit(cl.begin());cit!=cl.end();)
	      if (!((*cit)->Id()&cs[i]->m_id)) ++cit;
	      else cit=cl.erase(cit);
	    cl.push_back(ampl->Legs().back());
	    inc|=cs[i]->m_id;
	  }
	  SortFlavours(cl,s_usefmm?&fmm:NULL);
	  (*fit)->Delete();
	  fit=fl.erase(fit);
	  fl.insert(fit,cl.begin(),cl.end());
	  ampl->Legs().pop_back();
	  break;
	}
  }
  for (size_t i(0);i<ampl->NIn();++i) {
    il[i]->SetFlav(il[i]->Flav().Bar());
    ampl->Legs()[i]=il[i];
  }
  for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
    ampl->Legs()[i]=fl[i-ampl->NIn()];
}

std::string Process_Base::GenerateName(const Cluster_Amplitude *ampl)
{
  char nii[3], nfi[3];
  if (sprintf(nii,"%i",(int)ampl->NIn())<=0)
    THROW(fatal_error,"Conversion error");
  if (sprintf(nfi,"%i",(int)(ampl->Legs().size()-ampl->NIn()))<=0)
    THROW(fatal_error,"Conversion error");
  std::string name(nii+std::string("_")+nfi);
  for (size_t i(0);i<ampl->NIn();++i) 
    name+="__"+ampl->Leg(i)->Flav().Bar().IDName();
  DecayInfo_Vector decs(ampl->Decays());
  std::sort(decs.begin(),decs.end(),Order_NDecay());
  for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
    std::string op, cl;
    for (size_t j(0);j<decs.size();++j) {
      Int_Vector ids(ID(decs[j]->m_id));
      if (ampl->Leg(i)->Id()==(1<<ids.front()))
	op+=ToString(decs[j]->m_fl)+"[";
      else if (ampl->Leg(i)->Id()==(1<<ids.back())) cl+="]";
    }
    name+="__"+op+ampl->Leg(i)->Flav().IDName()+cl;
  }
  return name;
}

std::string Process_Base::GenerateName(const NLO_subevt *sub,const size_t &nin)
{
  char nii[3], nfi[3];
  if (sprintf(nii,"%i",(int)nin)<=0)
    THROW(fatal_error,"Conversion error");
  if (sprintf(nfi,"%i",(int)(sub->m_n-nin))<=0)
    THROW(fatal_error,"Conversion error");
  std::string name(nii+std::string("_")+nfi);
  for (size_t i(0);i<sub->m_n;++i) name+="__"+sub->p_fl[i].IDName();
  return name;
}

void Process_Base::SetGenerator(ME_Generator_Base *const gen) 
{ 
  p_gen=gen; 
}

void Process_Base::SetShower(PDF::Shower_Base *const ps)
{
  p_shower=ps; 
}

void Process_Base::SetVariationWeights(SHERPA::Variation_Weights *const vw)
{
  if (m_variationweightsowned) {
    delete p_variationweights;
    m_variationweightsowned = false;
  }
  p_variationweights=vw;
  if (p_int->PSHandler() != NULL) p_int->PSHandler()->SetVariationWeights(vw);
}

void Process_Base::SetOwnedVariationWeights(SHERPA::Variation_Weights *vw)
{
  SetVariationWeights(vw);
  if (vw) {
    m_variationweightsowned = true;
  }
}

void Process_Base::FillOnshellConditions()
{
  if (!Selector()) return;
  Subprocess_Info info(m_pinfo.m_ii);
  info.Add(m_pinfo.m_fi);
  for(size_t i=0;i<m_decins.size();i++)
    if (m_decins[i]->m_osd) Selector()->AddOnshellCondition
      (PSId(m_decins[i]->m_id),sqr(m_decins[i]->m_fl.Mass()));
}

void Process_Base::FillAmplitudes(std::vector<METOOLS::Spin_Amplitudes>& amp,
                                  std::vector<std::vector<Complex> >& cols)
{
  THROW(fatal_error, "Virtual function called.");
}

void Process_Base::SetSelector(const Selector_Key &key)
{
  if (IsMapped()) return;
  if (p_selector==NULL) p_selector = new Combined_Selector(p_int);
  p_selector->Initialize(key);
}

bool Process_Base::Trigger(const Vec4D_Vector &p)
{
  if (IsMapped() && LookUp()) return Selector()->Result();
  return Selector()->Trigger(p);
}
 
NLO_subevtlist *Process_Base::GetSubevtList()
{
  return NULL;
}

NLO_subevtlist *Process_Base::GetRSSubevtList()
{
  return NULL;
}

void Process_Base::InitCuts(Cut_Data *const cuts)
{
  cuts->Init(m_nin,m_flavs);
}

void Process_Base::BuildCuts(Cut_Data *const cuts)
{
  if (IsMapped() && LookUp()) return;
  Selector()->BuildCuts(cuts);
}

void Process_Base::SetRBMap(Cluster_Amplitude *ampl)
{
}

void Process_Base::InitPSHandler
(const double &maxerr,const std::string eobs,const std::string efunc)
{
  p_int->SetPSHandler(new Phase_Space_Handler(p_int,maxerr));
  p_int->PSHandler()->SetVariationWeights(p_variationweights);
  if (eobs!="") p_int->PSHandler()->SetEnhanceObservable(eobs);
  if (efunc!="") p_int->PSHandler()->SetEnhanceFunction(efunc);
} 

double Process_Base::LastPlus()
{
  if (IsGroup()) {
    double last=0.0;
    for (size_t i(0);i<Size();++i)
      last+=(*this)[i]->LastPlus();
    return last;
  }
  double last(Last());
  return last>0.0?last:0.0;
}

double Process_Base::LastMinus()
{
  if (IsGroup()) {
    double last=0.0;
    for (size_t i(0);i<Size();++i)
      last+=(*this)[i]->LastMinus();
    return last;
  }
  double last(Last());
  return last<0.0?last:0.0;
}

void Process_Base::FillProcessMap(NLOTypeStringProcessMap_Map *apmap)
{
  p_apmap=apmap;
  if (IsGroup()) {
    for (size_t i(0);i<Size();++i) (*this)[i]->FillProcessMap(apmap);
  }
  else {
    nlo_type::code nlot(m_pinfo.m_fi.m_nloqcdtype);
    std::string fname(m_name);
    size_t pos=fname.find("EW");
    if (pos!=std::string::npos) fname=fname.substr(0,pos-2);
    pos=fname.find("QCD");
    if (pos!=std::string::npos) fname=fname.substr(0,pos-2);
    if (nlot&nlo_type::vsub) nlot=nlo_type::vsub;
    if (nlot&nlo_type::rsub) nlot=nlo_type::rsub;
    if (apmap->find(nlot)==apmap->end())
      (*apmap)[nlot] = new StringProcess_Map();
    StringProcess_Map *cmap((*apmap)[nlot]);
    if (msg_LevelIsDebugging() && cmap->find(fname) != cmap->end()) {
      Process_Base* old = (*cmap)[fname];
      msg_Out()
        << METHOD << "(): replacing '" << m_name << "' "
        << Demangle(typeid(*old).name())
        << " -> "
        << Demangle(typeid(*this).name())
        << "\n";
    }
    (*cmap)[fname]=this;
  }
}

size_t Process_Base::SetMCMode(const size_t mcmode)
{
  size_t cmcmode(m_mcmode);
  m_mcmode=mcmode;
  if (IsGroup())
    for (size_t i(0);i<Size();++i)
      (*this)[i]->SetMCMode(mcmode);
  return cmcmode;
}

size_t Process_Base::SetClusterMode(const size_t cmode)
{
  size_t ccmode(m_cmode);
  m_cmode=cmode;
  if (IsGroup())
    for (size_t i(0);i<Size();++i)
      (*this)[i]->SetClusterMode(cmode);
  return ccmode;
}

void Process_Base::SetSProc(Process_Base *sproc)
{
  p_sproc=sproc;
  if (IsGroup())
    for (size_t i(0);i<Size();++i)
      (*this)[i]->SetSProc(sproc);
}

void Process_Base::SetBBarMC(BBar_Multi_Channel *mc)
{
  p_mc=mc;
  if (IsGroup())
    for (size_t i(0);i<Size();++i)
      (*this)[i]->SetBBarMC(mc);
}

int Process_Base::NaiveMapping(Process_Base *proc) const
{
  DEBUG_FUNC(Name()<<" -> "<<proc->Name());
  const Vertex_Table *vt(s_model->VertexTable());
  std::map<Flavour,Flavour> fmap;
  std::vector<Flavour> curf(m_flavs);
  for (size_t i(0);i<curf.size();++i) fmap[curf[i]]=proc->m_flavs[i];
  for (std::map<Flavour,Flavour>::const_iterator
	 fit(fmap.begin());fit!=fmap.end();++fit)
    DEBUG_VAR(fit->first<<" -> "<<fit->second);
  for (size_t i(0);i<curf.size();++i) {
    Flavour cf(curf[i]), mf(fmap[cf]);
    if (cf==mf) continue;
    const Vertex_List &vlc(vt->find(cf)->second);
    const Vertex_List &vlm(vt->find(mf)->second);
    DEBUG_VAR(cf<<" "<<vlc.size());
    DEBUG_VAR(mf<<" "<<vlm.size());
    if (vlc.size()!=vlm.size()) return 0;
    for (size_t j(0);j<vlc.size();++j) {
      msg_Indent();
      DEBUG_VAR(*vlc[j]);
      bool match(false);
      for (size_t k(0);k<vlm.size();++k) {
	msg_Indent();
	DEBUG_VAR(*vlm[k]);
	if (vlm[k]->Compare(vlc[j])==0) {
	  msg_Indent();
	  for (int n=1;n<vlc[j]->NLegs();++n) {
	    std::map<Flavour,Flavour>::
	      const_iterator cit(fmap.find(vlc[j]->in[n]));
	    if (cit==fmap.end()) {
	      DEBUG_VAR(vlc[j]->in[n]<<" -> "<<vlm[k]->in[n]);
	      if (vlc[j]->in[n].Mass()!=vlm[k]->in[n].Mass() ||
		  vlc[j]->in[n].Width()!=vlm[k]->in[n].Width()) {
		msg_Debugging()<<"m_"<<vlc[j]->in[n]
			       <<" = "<<vlc[j]->in[n].Mass()
			       <<" / m_"<<vlm[k]->in[n]
			       <<" = "<<vlm[k]->in[n].Mass()<<"\n";
		msg_Debugging()<<"w_"<<vlc[j]->in[n]
			       <<" = "<<vlc[j]->in[n].Width()
			       <<" / w_"<<vlm[k]->in[n]
			       <<" = "<<vlm[k]->in[n].Width()<<"\n";
		return 0;
	      }
	      fmap[vlc[j]->in[n]]=vlm[k]->in[n];
	      curf.push_back(vlc[j]->in[n]);
	    }
	    else if (cit->second!=vlm[k]->in[n]) {
	      DEBUG_VAR(cit->second<<" "<<vlm[k]->in[n]);
	      return 0;
	    }
	  }
	  DEBUG_VAR(*vlc[j]);
	  DEBUG_VAR(*vlm[k]);
	  match=true;
	  break;
	}
      }
      if (!match) return 0;
    }
  }
  DEBUG_VAR("OK");
  return 1;
}

std::string Process_Base::ShellName(std::string name) const
{
  if (name.length()==0) name=m_name;
  for (size_t i(0);(i=name.find('-',i))!=std::string::npos;name.replace(i,1,"m"));
  for (size_t i(0);(i=name.find('+',i))!=std::string::npos;name.replace(i,1,"p"));
  for (size_t i(0);(i=name.find('~',i))!=std::string::npos;name.replace(i,1,"x"));
  for (size_t i(0);(i=name.find('(',i))!=std::string::npos;name.replace(i,1,"_"));
  for (size_t i(0);(i=name.find(')',i))!=std::string::npos;name.replace(i,1,"_"));
  for (size_t i(0);(i=name.find('[',i))!=std::string::npos;name.replace(i,1,"I"));
  for (size_t i(0);(i=name.find(']',i))!=std::string::npos;name.replace(i,1,"I"));
  return name;
}
