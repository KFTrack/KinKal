#ifndef KinKal_HitInfo_hh
#define KinKal_HitInfo_hh
#include <vector>
namespace KinKal {
  struct HitInfo {
    enum htype{straw=0,scint=1,constraint=2,unknown};
    HitInfo(){};
    ~HitInfo(){};
    Int_t active_, type_, ambig_;
    Float_t time_, resid_, residvar_, fitchi_;
    static std::string leafnames() { return std::string("active/i:type/i:ambig/i:time/f:resid/f:residvar/f:fitchi/f"); }
  };
  typedef std::vector<HitInfo> KKHIV;
}
#endif
