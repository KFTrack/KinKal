#ifndef KinKal_HitInfo_hh
#define KinKal_HitInfo_hh
#include <vector>
namespace KinKal {
  struct HitInfo {
    HitInfo(){};
    ~HitInfo(){};
    Int_t active_;
    Float_t time_, resid_, residvar_, fitchi_;
    static std::string leafnames() { return std::string("active/i:time/f:resid/f:residvar/f:fitchi/f"); }
  };
  typedef std::vector<HitInfo> KKHIV;
}
#endif
