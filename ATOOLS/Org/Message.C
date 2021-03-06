#include "ATOOLS/Org/Message.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/My_MPI.H"

#include <sys/stat.h>
#include <iterator>

namespace ATOOLS {
  Message *msg(NULL);
}

using namespace ATOOLS;

std::ostream &ATOOLS::operator<<(std::ostream &str,const bm::code modifier) 
{
  switch (modifier) {
  case bm::back: return msg->Modifiable()?str<<"\b":str<<" \\b ";
  case bm::cr:   return msg->Modifiable()?str<<"\r":str<<"\n";
  case bm::bell: return msg->Modifiable()?str<<"\a":str<<" \\a ";
  case bm::none: return str;
  }
  return str;
}

std::ostream &ATOOLS::operator<<(std::ostream &str,const om::code modifier) 
{
  if (!msg->Modifiable()) return str;
  switch (modifier) {
#ifdef USING__COLOUR
  case om::reset:    return str<<"\033[0m";
  case om::bold:     return str<<"\033[1m";
  case om::underln:  return str<<"\033[4m";
  case om::blink:    return str<<"\033[5m";
  case om::blackbg:  return str<<"\033[7m";
  case om::red:      return str<<"\033[31m";
  case om::green:    return str<<"\033[32m";
  case om::brown:    return str<<"\033[33m";
  case om::blue:     return str<<"\033[34m";
  case om::violet:   return str<<"\033[35m";
  case om::lblue:    return str<<"\033[36m";
  case om::grey:     return str<<"\033[37m";
  case om::redbg:    return str<<"\033[41m";
  case om::greenbg:  return str<<"\033[42m";
  case om::brownbg:  return str<<"\033[43m";
  case om::bluebg:   return str<<"\033[44m";
  case om::violetbg: return str<<"\033[45m";
  case om::lbluebg:  return str<<"\033[46m";
  case om::greybg:   return str<<"\033[47m";
  case om::none:     return str;
#else
  default: return str;
#endif
  }
  return str;
}
 
std::ostream &ATOOLS::operator<<(std::ostream &str,const mm modifier)
{
  if (!msg->Modifiable()) return str;
  switch (modifier.m_code) {
#ifdef USING__COLOUR
  case mm::up:    return str<<"\033["<<modifier.m_num<<"A";
  case mm::down:  return str<<"\033["<<modifier.m_num<<"B";
  case mm::right: return str<<"\033["<<modifier.m_num<<"C";
  case mm::left:  return str<<"\033["<<modifier.m_num<<"D";
  case mm::none:  return str;
#else
  default: return str;
#endif
  }
  return str;
}

std::ostream &ATOOLS::operator<<(std::ostream &str,const tm::code modifier) 
{
  if (!msg->Modifiable()) return str;
  switch (modifier) {
#ifdef USING__COLOUR
  case tm::curon:  return str<<"\033[?25h";
  case tm::curoff: return str<<"\033[?25l";
  case tm::none:   return str;
#else
  default: return str;
#endif
  }
  return str;
}


indentbuf::indentbuf(std::streambuf* basebuf) :
  m_basebuf(basebuf), m_indent(0), at_start(true)
{
}

indentbuf::~indentbuf()
{
}

void indentbuf::Indent(size_t i)
{
  m_indent+=i;
}

void indentbuf::DeIndent(size_t i)
{
  if (m_indent>=i) m_indent-=i;
}

std::streambuf::int_type indentbuf::overflow(int_type ch)
{
  if (traits_type::eq_int_type(ch, traits_type::to_int_type('\r'))) {
  }
  if (ch == traits_type::eof())
    return traits_type::not_eof(ch);

  if (traits_type::not_eof(ch)) {
    if (at_start)
      for (size_t i = 0; i < m_indent; ++i)
        m_basebuf->sputc(traits_type::to_char_type(' '));
    m_basebuf->sputc(traits_type::to_char_type(ch));
    if (traits_type::eq_int_type(ch, traits_type::to_int_type('\n')))
      at_start = true;
    else
      at_start = false;
  }
  return ch; 
}

Message::Message() : m_buf(std::cout.rdbuf()), p_log(NULL)
{      
  p_output = &std::cout;
  p_error = &std::cerr;
  p_no = new std::ofstream("/dev/null",std::ios::app);
  m_logfile = "";
  m_file = 0;
  m_level = 0;
  m_modifiable = true;
  m_mpimode = 0;
}

Message::~Message() 
{      
  p_output->rdbuf(m_buf.BaseBuf());
  if (p_log) delete p_log;
  delete p_no;
}

void Message::Init(const std::string& level,const std::string &logfile) 
{ 
  if (logfile!="") {
    m_logfile = logfile;
    p_log = new std::ofstream(logfile.c_str(),std::ios::app);
    p_output=p_log;
  }
  m_buf.SetBaseBuf(p_output->rdbuf());
  p_output->rdbuf(&m_buf);

  Data_Reader dr("|","[","!");
  dr.AddLineSeparator("]");
  dr.AddComment("#");
  dr.SetString(level);
  std::vector<std::vector<std::string> > svv;
  dr.MatrixFromString(svv);

  m_level=2;
  for (size_t i=0;i<svv.size();++i) {
    if (svv[i].size()==1) m_level=ATOOLS::ToType<int>(svv[i][0]);
    else if (svv[i].size()==2) {
      int lev=ToType<int>(svv[i][1]);
      if (lev&1) m_contextevents.insert(svv[i][0]);
      if (lev&2) m_contextinfo.insert(svv[i][0]);
      if (lev&4) m_contexttracking.insert(svv[i][0]);
      if (lev&8) m_contextdebugging.insert(svv[i][0]);
      if (lev&32) m_contextiodebugging.insert(svv[i][0]);
    }
    else {
      Out()<<METHOD<<": Do not understand output specification: "<<std::endl;
      copy(svv[i].begin(), svv[i].end(),
           std::ostream_iterator<std::string>(Out(), " "));
    }
  }
}

void Message::SetStandard() 
{
  if (m_file==1) delete p_output;
  p_output = &std::cout;
  p_error = &std::cerr;
  m_file = 0;
}

void Message::SetPrecision(const int precision) 
{
  if (p_output) p_output->precision(precision);
}

std::ostream &Message::Out() const 
{ 
#ifdef USING__MPI
  if (!m_mpimode && mpi->Rank()) return *p_no;
#endif
  return *p_output; 
}

std::ostream &Message::Error() const     
{ 
#ifdef USING__MPI
  if (!m_mpimode && 
      mpi->Rank()) return *p_no;
#endif
  if (m_level >= 0) return *p_output; 
  return *p_no; 
}

std::ostream &Message::Events() const    
{ 
#ifdef USING__MPI
  if (!m_mpimode && 
      mpi->Rank()) return *p_no;
#endif
  if (m_level & 1) return *p_output; 
  return *p_no;  
}

std::ostream &Message::Info() const      
{ 
#ifdef USING__MPI
  if (!m_mpimode && 
      mpi->Rank()) return *p_no;
#endif
  if (m_level & 2) return *p_output; 
  return *p_no;  
}

std::ostream &Message::Tracking() const  
{ 
#ifdef USING__MPI
  if (!m_mpimode && 
      mpi->Rank()) return *p_no;
#endif
  if (m_level & 4) return *p_output; 
  return *p_no;  
}

std::ostream &Message::Debugging() const 
{ 
#ifdef USING__MPI
  if (!m_mpimode && mpi->Rank()) return *p_no;
#endif
  if (m_level & 8) return *p_output; 
  return *p_no;  
}

std::ostream &Message::IODebugging() const
{
#ifdef USING__MPI
  if (!m_mpimode && mpi->Rank()) return *p_no;
#endif
  if (m_level & 32) return *p_output;
  return *p_no;
}

std::string Message::ExtractMethodName(std::string cmethod) const   
{ 
  for (size_t pos(cmethod.find(", "));
       pos!=std::string::npos;pos=cmethod.find(", ")) cmethod.erase(pos+1,1);
  for (size_t pos(cmethod.find("> >"));
       pos!=std::string::npos;pos=cmethod.find("> >")) cmethod.erase(pos+1,1);
  std::string cclass("<no class>"), method("<no method>");
  cmethod=cmethod.substr(0,ATOOLS::Min(cmethod.length(),cmethod.find("(")));
  size_t pos;
  while ((pos=cmethod.find(" "))!=std::string::npos) 
    cmethod=cmethod.substr(pos+1);
  pos=cmethod.find("::");
  for (size_t bpos(cmethod.find("<"));pos!=std::string::npos && pos<bpos;bpos-=pos+2) {
    cclass=cmethod.substr(0,pos);
    cmethod=cmethod.substr(pos+2);
    pos=cmethod.rfind("::");
    method=cmethod.substr(0,ATOOLS::Min(cmethod.length(),pos));
  }
  if (cclass=="<no class>") return cmethod;
  return cclass+"::"+cmethod;
}

bool Message::LevelIsEvents(const std::string& context) const
{
  for (std::set<std::string>::reverse_iterator rit=m_contextevents.rbegin();
       rit!=m_contextevents.rend(); ++rit) {
    if (context.find(*rit)!=std::string::npos) return true;
  }
  return false;
}

bool Message::LevelIsInfo(const std::string& context) const
{
  for (std::set<std::string>::reverse_iterator rit=m_contextinfo.rbegin();
       rit!=m_contextinfo.rend(); ++rit) {
    if (context.find(*rit)!=std::string::npos) return true;
  }
  return false;
}

bool Message::LevelIsTracking(const std::string& context) const
{
  for (std::set<std::string>::reverse_iterator rit=m_contexttracking.rbegin();
       rit!=m_contexttracking.rend(); ++rit) {
    if (context.find(*rit)!=std::string::npos) return true;
  }
  return false;
}

bool Message::LevelIsDebugging(const std::string& context) const
{
  for (std::set<std::string>::reverse_iterator rit=m_contextdebugging.rbegin();
       rit!=m_contextdebugging.rend(); ++rit) {
    if (context.find(*rit)!=std::string::npos) return true;
  }
  return false;
}

bool Message::LevelIsIODebugging(const std::string& context) const
{
  for (std::set<std::string>::reverse_iterator rit=m_contextiodebugging.rbegin();
       rit!=m_contextiodebugging.rend(); ++rit) {
    if (context.find(*rit)!=std::string::npos) return true;
  }
  return false;
}
