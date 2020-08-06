#ifndef KinKal_KKHitInfo_hh
#define KinKal_KKHitInfo_hh
#include <vector>
namespace KinKal {
  struct KKHitInfo {
    KKHitInfo(){};
    ~KKHitInfo(){};
    Int_t active_;
    Float_t time_, resid_, residvar_, fitchi_;
    static std::string leafnames() { return std::string("active/i:time/f:resid/f:residvar/f:fitchi/f"); }
  };
  typedef std::vector<KKHitInfo> KKHIV;
}
#endif
